/*
 * HAL format 3: changing an alignment's genomes.  See hal3Update.h.
 */
#include "hal3Update.h"
#include "liftmap.h"
#include "twobit64.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <functional>
#include <set>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

using namespace std;

namespace hal3 {

namespace {

typedef vector<pair<string, int64_t>> SeqTable;

struct LmapGuard {
    lmap_file *f = nullptr;
    ~LmapGuard() {
        if (f) lmap_close(f);
    }
};

lmap_file *openHgn(const Store &store, const string &member) {
    if (!store.has(member)) throw runtime_error(store.path() + " has no " + member);
    lmap_file *f = lmap_open_io(store.openMember(member), 1);
    if (!f) throw runtime_error("cannot read " + member + " of " + store.path());
    return f;
}

SeqTable table(lmap_file *f, int axis) {
    SeqTable t;
    for (uint32_t i = 0; i < lmap_n_members(f, axis); i++) {
        const lmap_member *m = lmap_member_at(f, axis, i);
        t.push_back({m->name, (int64_t)m->length});
    }
    return t;
}

string meta(lmap_file *f, const char *key) {
    const char *v = lmap_meta_get(f, key);
    return v ? v : "";
}

/* Does the genome file f fit genome `name` with parent `parent` ("" for the root)? */
bool fits(lmap_file *f, const string &name, const string &parent) {
    if (meta(f, "hal.genome") != name || meta(f, "hal.parent") != parent) return false;
    lmap_file_stats st;
    lmap_get_stats(f, &st);
    bool branch = st.runs > 0 || lmap_n_members(f, 1) > 0;
    return parent.empty() ? !branch : true;
}

/* A copy of genome file src for genome `name` with parent `parent`: the names in its
 * metadata changed, and for the root ("" parent) its branch dropped. */
void rewriteHgn(lmap_file *src, const string &name, const string &parent, const string &out) {
    lmap_builder *b = lmap_builder_open(out.c_str(), HGN_PROFILE, 0, nullptr);
    if (!b) throw runtime_error("cannot write " + out);
    auto ok = [&](int rc) {
        if (rc < 0) {
            string e = lmap_builder_error(b);
            lmap_builder_abort(b);
            throw runtime_error(out + ": " + e);
        }
    };
    ok(lmap_builder_set_params(b, 8192, 1000000, LMAP_CODEC_DEFLATE));
    ok(lmap_builder_set_order(b, LMAP_ORDER_A));
    ok(lmap_builder_set_presorted(b, 1));
    for (uint32_t i = 0; i < lmap_meta_count(src); i++) {
        const char *k, *v;
        if (lmap_meta_at(src, i, &k, &v) != 0) continue;
        string key = k;
        if (key == "hal.genome" || key == "hal.parent" || key == "axis.a" || key == "axis.b") continue;
        ok(lmap_builder_set_meta(b, k, v));
    }
    ok(lmap_builder_set_meta(b, "axis.a", ("genome " + name).c_str()));
    ok(lmap_builder_set_meta(b, "hal.genome", name.c_str()));
    if (!parent.empty()) {
        ok(lmap_builder_set_meta(b, "axis.b", ("genome " + parent).c_str()));
        ok(lmap_builder_set_meta(b, "hal.parent", parent.c_str()));
    }
    for (uint32_t i = 0; i < lmap_n_members(src, 0); i++) {
        const lmap_member *m = lmap_member_at(src, 0, i);
        ok(lmap_builder_add_member(b, 0, m->name, m->length));
    }
    if (!parent.empty()) {
        for (uint32_t i = 0; i < lmap_n_members(src, 1); i++) {
            const lmap_member *m = lmap_member_at(src, 1, i);
            ok(lmap_builder_add_member(b, 1, m->name, m->length));
        }
        for (uint32_t i = 0; i < lmap_n_members(src, 0); i++) {
            lmap_cursor *c = lmap_cursor_open(src, 0, i, nullptr, 0, 0);
            if (!c) ok(-1);
            lmap_hit h;
            int rc = lmap_cursor_seek(c, 0, (int64_t)lmap_member_at(src, 0, i)->length);
            while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
                rc = lmap_builder_add_run(b, h.a_member, h.b_member, h.a, h.b, h.len, h.strand);
                if (rc < 0) break;
                rc = 0;
            }
            lmap_cursor_close(c);
            ok(rc);
        }
    }
    lmap_build_stats st;
    if (lmap_builder_close(b, &st) != 0) throw runtime_error(out + ": " + st.error);
}

/* Copy an lmap_io to a file; returns the CRC-32 and sets *bytes. */
uint32_t copyToFile(lmap_io *io, const string &path, uint64_t *bytes) {
    FILE *f = fopen(path.c_str(), "wb");
    if (!f) throw runtime_error("cannot write " + path + ": " + strerror(errno));
    int64_t n = lmap_io_size(io);
    vector<unsigned char> buf((size_t)4 << 20);
    uLong c = crc32(0L, Z_NULL, 0);
    bool good = true;
    for (int64_t at = 0; at < n && good; at += (int64_t)buf.size()) {
        int64_t k = min<int64_t>((int64_t)buf.size(), n - at);
        good = lmap_pread(io, buf.data(), at, k) == 0 && fwrite(buf.data(), 1, (size_t)k, f) == (size_t)k;
        c = crc32(c, buf.data(), (uInt)k);
    }
    good = (fclose(f) == 0) && good;
    if (!good) throw runtime_error("cannot write " + path);
    *bytes = (uint64_t)n;
    return (uint32_t)c;
}

string quoteNewick(const string &name) {
    if (name.find_first_of("()[]',:; \t") == string::npos) return name;
    string q = "'";
    for (char c : name) q += c == '\'' ? string("''") : string(1, c);
    return q + "'";
}

string fmtLength(double d) {
    char buf[64];
    for (int prec = 6; prec <= 17; prec++) {
        snprintf(buf, sizeof buf, "%.*g", prec, d);
        if (strtod(buf, nullptr) == d) break;
    }
    return buf;
}

/* A scratch directory beside path, removed with its contents. */
struct ScratchDir {
    string dir;
    explicit ScratchDir(const string &path) : dir(path + ".update-tmp") {
        if (mkdir(dir.c_str(), 0777) != 0)
            throw runtime_error("cannot create " + dir + (errno == EEXIST ? " (left by an interrupted update? remove it)" : ""));
    }
    ~ScratchDir() {
        for (const string &f : files) unlink(f.c_str());
        rmdir(dir.c_str());
    }
    string file(const string &stem) {
        files.push_back(dir + "/" + stem);
        return files.back();
    }
    vector<string> files;
};
}

/* ---------------------------------------------------------------- the model */

Update::Update(const string &path) : _path(path) {
    _in = Store::open(path);
    _manifest = Manifest::parse(_in->read(MANIFEST_NAME), path);
    for (const GenomeEntry &e : _manifest.genomes) {
        Node n;
        n.name = e.name;
        n.parent = e.parent;
        n.branchLength = e.branchLength;
        n.hgn = Source{_in.get(), e.hgnPath};
        n.dna = Source{_in.get(), e.dnaPath};
        _nodes[e.name] = n;
        if (e.parent.empty()) _root = e.name;
    }
    for (const GenomeEntry &e : _manifest.genomes)
        if (!e.parent.empty()) need(e.parent).children.push_back(e.name);
}

Update::~Update() {}

Update::Node &Update::need(const string &name) {
    auto it = _nodes.find(name);
    if (it == _nodes.end()) throw runtime_error("no genome " + name + " in the alignment being updated");
    return it->second;
}

void Update::detach(const string &name) {
    Node &n = need(name);
    auto p = _nodes.find(n.parent);
    if (p != _nodes.end()) {
        auto &ch = p->second.children;
        ch.erase(std::remove(ch.begin(), ch.end(), name), ch.end());
    }
}

void Update::dropSubtree(const string &name) {
    vector<string> children = need(name).children;
    for (const string &c : children) dropSubtree(c);
    _nodes.erase(name);
}

void Update::remove(const vector<string> &genomes) {
    for (const string &g : genomes) {
        detach(g);
        _nodes.erase(g);
        if (g == _root) _root.clear();
        _treeChanged = true;
        _log.push_back("removed " + g);
    }
}

void Update::prune(const string &genome) {
    detach(genome);
    size_t before = _nodes.size();
    dropSubtree(genome);
    if (_nodes.find(_root) == _nodes.end()) _root.clear();
    _treeChanged = true;
    _log.push_back("pruned " + genome + " (" + to_string(before - _nodes.size()) + " genomes)");
}

void Update::extract(const string &genome) {
    Node &g = need(genome);
    set<string> keepSet;
    vector<string> stack{genome};
    while (!stack.empty()) {
        string n = stack.back();
        stack.pop_back();
        keepSet.insert(n);
        for (const string &c : need(n).children) stack.push_back(c);
    }
    size_t before = _nodes.size();
    for (auto it = _nodes.begin(); it != _nodes.end();) it = keepSet.count(it->first) ? next(it) : _nodes.erase(it);
    g.parent.clear();
    g.branchLength = -1;
    _root = genome;
    _treeChanged = true;
    _log.push_back("extracted the subtree of " + genome + " (" + to_string(_nodes.size()) + " genomes kept, " +
                   to_string(before - _nodes.size()) + " removed)");
}

void Update::keep(const vector<string> &genomes) {
    if (genomes.empty()) throw runtime_error("--keep needs genomes");
    auto path = [&](const string &g) {   /* g, its parent, ... up to the root */
        vector<string> p;
        for (string x = g; !x.empty(); x = need(x).parent) {
            p.push_back(x);
            if (_nodes.find(need(x).parent) == _nodes.end()) break;
        }
        return p;
    };
    vector<string> first = path(genomes[0]);
    string lca = genomes[0];
    for (const string &g : genomes) {
        vector<string> p = path(g);
        set<string> ps(p.begin(), p.end());
        auto it = find_if(first.begin(), first.end(), [&](const string &x) { return ps.count(x) > 0; });
        if (it == first.end()) throw runtime_error("--keep: " + g + " is not in the same tree as " + genomes[0]);
        if (find(first.begin(), first.end(), lca) < it) lca = *it;
    }
    set<string> keepSet;
    for (const string &g : genomes)
        for (const string &x : path(g)) {
            keepSet.insert(x);
            if (x == lca) break;
        }
    size_t before = _nodes.size();
    for (auto it = _nodes.begin(); it != _nodes.end();) it = keepSet.count(it->first) ? next(it) : _nodes.erase(it);
    for (auto &n : _nodes) {
        auto &ch = n.second.children;
        ch.erase(std::remove_if(ch.begin(), ch.end(), [&](const string &c) { return !keepSet.count(c); }), ch.end());
    }
    Node &root = need(lca);
    root.parent.clear();
    root.branchLength = -1;
    _root = lca;
    _treeChanged = true;
    _log.push_back("kept " + to_string(_nodes.size()) + " genomes (the spanning tree of " + to_string(genomes.size()) +
                   ", rooted at " + lca + "), removed " + to_string(before - _nodes.size()));
}

void Update::graft(const string &donorPath) {
    _donors.push_back(Store::open(donorPath));
    const Store *d = _donors.back().get();
    Manifest m = Manifest::parse(d->read(MANIFEST_NAME), donorPath);
    if (m.genomes.empty()) throw runtime_error(donorPath + " has no genomes");
    const string &r = m.genomes[0].name;
    need(r);
    int added = 0;
    for (size_t i = 1; i < m.genomes.size(); i++) {
        const GenomeEntry &e = m.genomes[i];
        if (_nodes.count(e.name)) throw runtime_error("--graft: genome " + e.name + " of " + donorPath + " is already here");
        Node n;
        n.name = e.name;
        n.parent = e.parent;
        n.branchLength = e.branchLength;
        n.hgn = Source{d, e.hgnPath};
        n.dna = Source{d, e.dnaPath};
        _nodes[e.name] = n;
        need(e.parent).children.push_back(e.name);
        added++;
    }
    _treeChanged = true;
    _log.push_back("grafted " + to_string(added) + " genomes of " + donorPath + " below " + r);
}

void Update::replace(const string &donorPath, const vector<string> &genomes) {
    _donors.push_back(Store::open(donorPath));
    const Store *d = _donors.back().get();
    Manifest m = Manifest::parse(d->read(MANIFEST_NAME), donorPath);
    for (const string &g : genomes) {
        Node &n = need(g);
        const GenomeEntry *e = m.find(g);
        if (!e) throw runtime_error("--replace: " + donorPath + " has no genome " + g);
        if (e->parent.empty() && !n.parent.empty())
            throw runtime_error("--replace: " + g + " is the root of " + donorPath + ", so its file has no branch to put it below " +
                                n.parent);
        n.hgn = Source{d, e->hgnPath};
        n.dna = Source{d, e->dnaPath};
        if (e->parent != n.parent) {
            /* the donor places it elsewhere: move it there */
            detach(g);
            n.parent = e->parent;
            auto p = _nodes.find(e->parent);
            if (p != _nodes.end()) p->second.children.push_back(g);
            _treeChanged = true;
        }
        if (e->branchLength != n.branchLength) _treeChanged = true;
        n.branchLength = e->branchLength;
        _log.push_back("replaced " + g + " from " + donorPath);
    }
}

void Update::rename(const string &from, const string &to) {
    if (_nodes.count(to)) throw runtime_error("--rename: there is already a genome " + to);
    Node n = need(from);
    detach(from);
    _nodes.erase(from);
    n.name = to;
    auto p = _nodes.find(n.parent);
    if (p != _nodes.end()) p->second.children.push_back(to);
    for (const string &c : n.children)
        if (_nodes.count(c)) _nodes[c].parent = to;
    _nodes[to] = n;
    if (_root == from) _root = to;
    _treeChanged = true;
    _log.push_back("renamed " + from + " to " + to);
}

void Update::setBranchLength(const string &genome, double length) {
    Node &n = need(genome);
    if (n.parent.empty()) throw runtime_error("--branchLength: " + genome + " is the root");
    n.branchLength = length;
    _treeChanged = true;
    _log.push_back("branch length of " + genome + " set to " + fmtLength(length));
}

vector<string> Update::preorder() const {
    vector<string> out, stack;
    if (!_root.empty()) stack.push_back(_root);
    while (!stack.empty()) {
        string g = stack.back();
        stack.pop_back();
        out.push_back(g);
        const Node &n = _nodes.at(g);
        for (auto c = n.children.rbegin(); c != n.children.rend(); ++c)
            if (_nodes.count(*c)) stack.push_back(*c);
    }
    return out;
}

string Update::newick() const {
    function<string(const string &)> rec = [&](const string &g) -> string {
        const Node &n = _nodes.at(g);
        string s;
        vector<string> kids;
        for (const string &c : n.children)
            if (_nodes.count(c)) kids.push_back(rec(c));
        if (!kids.empty()) {
            s = "(";
            for (size_t i = 0; i < kids.size(); i++) s += (i ? "," : "") + kids[i];
            s += ")";
        }
        s += quoteNewick(g);
        if (!n.parent.empty() && n.branchLength >= 0) s += ":" + fmtLength(n.branchLength);
        return s;
    };
    return rec(_root) + ";";
}

/* ---------------------------------------------------------------- checking */

vector<string> Update::check() {
    vector<string> problems;
    vector<string> roots;
    for (auto &kv : _nodes)
        if (kv.second.parent.empty()) roots.push_back(kv.first);
    if (roots.size() != 1) {
        string r;
        for (const string &x : roots) r += " " + x;
        problems.push_back(roots.empty() ? "no root" : "more than one root:" + r);
        return problems;
    }
    _root = roots[0];
    vector<string> order = preorder();
    if (order.size() != _nodes.size()) {
        for (auto &kv : _nodes)
            if (find(order.begin(), order.end(), kv.first) == order.end())
                problems.push_back(kv.first + ": its parent " + kv.second.parent + " is not in the result");
        return problems;
    }
    /* each genome's parent table against its parent's own table */
    map<string, SeqTable> own;
    for (const string &g : order) {
        const Node &n = _nodes.at(g);
        LmapGuard f;
        try {
            f.f = openHgn(*n.hgn.store, n.hgn.member);
        } catch (exception &e) {
            problems.push_back(g + ": " + e.what());
            continue;
        }
        own[g] = table(f.f, 0);
        if (n.parent.empty()) continue;
        if (!own.count(n.parent)) continue;   /* the parent's own file failed: reported already */
        if (table(f.f, 1) != own[n.parent]) {
            if (lmap_n_members(f.f, 1) == 0)
                problems.push_back(g + ": its genome file has no branch (it was a root), so it cannot go below " + n.parent);
            else
                problems.push_back(g + ": its genome file maps it to a parent whose sequences are not " + n.parent + "'s");
        }
    }
    /* DNA against the sequence table */
    for (const string &g : order) {
        const Node &n = _nodes.at(g);
        if (n.dna.member.empty() || !own.count(g)) continue;
        if (!n.dna.store->has(n.dna.member)) {
            problems.push_back(g + ": " + n.dna.member + " is missing");
            continue;
        }
        lmap_io *io = n.dna.store->openMember(n.dna.member);
        tb_file *t = tb_open_io(io, 1);
        if (!t) {
            lmap_io_close(io);
            problems.push_back(g + ": cannot read " + n.dna.member);
            continue;
        }
        const SeqTable &st = own[g];
        bool same = tb_n_seqs(t) == st.size();
        for (uint32_t i = 0; same && i < st.size(); i++) same = st[i].first == tb_name(t, i) && tb_length(t, i) == st[i].second;
        tb_close(t);
        if (!same) problems.push_back(g + ": its DNA does not match its sequence table");
    }
    return problems;
}

/* ---------------------------------------------------------------- writing */

Manifest Update::finalManifest() const {
    Manifest m;
    m.newick = _treeChanged ? newick() : _manifest.newick;
    m.meta = _manifest.meta;
    for (const string &g : preorder()) {
        const Node &n = _nodes.at(g);
        GenomeEntry e;
        e.name = g;
        e.parent = n.parent;
        e.branchLength = n.parent.empty() ? -1 : n.branchLength;
        e.hgnPath = "genomes/" + genomeFileStem(g) + ".hgn";
        if (!n.dna.member.empty()) e.dnaPath = "genomes/" + genomeFileStem(g) + ".2bit";
        m.genomes.push_back(e);
    }
    return m;
}

namespace {
string problemsText(const vector<string> &p) {
    string s = "the result would not be a valid alignment:";
    for (const string &x : p) s += "\n  " + x;
    return s;
}
}

void Update::write(const string &out, bool directory) {
    vector<string> p = check();
    if (!p.empty()) throw runtime_error(problemsText(p));
    struct stat st;
    if (stat(out.c_str(), &st) == 0) throw runtime_error(out + " already exists");
    Manifest m = finalManifest();
    ScratchDir scratch(out);
    unique_ptr<PackWriter> w;
    if (directory) {
        for (const string &d : {out, out + "/genomes"})
            if (mkdir(d.c_str(), 0777) != 0) throw runtime_error("cannot create " + d + ": " + strerror(errno));
    } else {
        w = PackWriter::create(out);
    }
    auto put = [&](const string &member, lmap_io *io, uint64_t *bytes) -> uint32_t {
        if (w) {
            *bytes = (uint64_t)lmap_io_size(io);
            return w->putIo(member, io);
        }
        return copyToFile(io, out + "/" + member, bytes);
    };
    for (GenomeEntry &e : m.genomes) {
        const Node &n = _nodes.at(e.name);
        LmapGuard f;
        f.f = openHgn(*n.hgn.store, n.hgn.member);
        SeqTable t = table(f.f, 0);
        e.numSequences = t.size();
        e.length = 0;
        for (auto &s : t) e.length += (uint64_t)s.second;
        lmap_io *io;
        if (fits(f.f, e.name, e.parent)) {
            io = n.hgn.store->openMember(n.hgn.member);
        } else {
            string tmp = scratch.file(genomeFileStem(e.name) + ".hgn");
            rewriteHgn(f.f, e.name, e.parent, tmp);
            io = lmap_io_open_file(tmp.c_str());
            _log.push_back("rewrote the genome file of " + e.name);
        }
        e.hgnCrc = put(e.hgnPath, io, &e.hgnBytes);
        lmap_io_close(io);
        if (!e.dnaPath.empty()) {
            io = n.dna.store->openMember(n.dna.member);
            e.dnaCrc = put(e.dnaPath, io, &e.dnaBytes);
            lmap_io_close(io);
        }
    }
    if (w) {
        w->putBytes(MANIFEST_NAME, m.text());
        w->commit();
    } else {
        m.write(out + "/" + MANIFEST_NAME);
    }
}

void Update::writeInPlace() {
    vector<string> p = check();
    if (!p.empty()) throw runtime_error(problemsText(p));
    Manifest m = finalManifest();
    ScratchDir scratch(_path);
    /* what each member was in the input, to leave unchanged members alone */
    map<string, const GenomeEntry *> before;
    for (const GenomeEntry &e : _manifest.genomes) {
        before[e.hgnPath] = &e;
        if (!e.dnaPath.empty()) before[e.dnaPath] = &e;
    }
    unique_ptr<PackWriter> w;
    if (_in->isPack()) w = PackWriter::update(_path);
    vector<pair<string, string>> renames;   /* directory: temp file -> member path */
    auto put = [&](const string &member, lmap_io *io, uint64_t *bytes) -> uint32_t {
        if (w) {
            *bytes = (uint64_t)lmap_io_size(io);
            return w->putIo(member, io);
        }
        string tmp = scratch.file(to_string(renames.size()) + ".new");
        uint32_t c = copyToFile(io, tmp, bytes);
        renames.push_back({tmp, _path + "/" + member});
        return c;
    };
    set<string> finalMembers{MANIFEST_NAME};
    for (GenomeEntry &e : m.genomes) {
        const Node &n = _nodes.at(e.name);
        finalMembers.insert(e.hgnPath);
        if (!e.dnaPath.empty()) finalMembers.insert(e.dnaPath);
        LmapGuard f;
        f.f = openHgn(*n.hgn.store, n.hgn.member);
        SeqTable t = table(f.f, 0);
        e.numSequences = t.size();
        e.length = 0;
        for (auto &s : t) e.length += (uint64_t)s.second;
        bool fit = fits(f.f, e.name, e.parent);
        /* unchanged: the input's own member, under its own name, as it is */
        auto was = before.find(e.hgnPath);
        auto from = before.find(n.hgn.member);
        if (fit && n.hgn.store == _in.get() && n.hgn.member == e.hgnPath && was != before.end()) {
            e.hgnBytes = was->second->hgnBytes;
            e.hgnCrc = was->second->hgnCrc;
        } else if (w && fit && n.hgn.store == _in.get() && from != before.end()) {
            w->link(n.hgn.member, e.hgnPath);
            e.hgnBytes = from->second->hgnBytes;
            e.hgnCrc = from->second->hgnCrc;
        } else {
            lmap_io *io;
            if (fit) {
                io = n.hgn.store->openMember(n.hgn.member);
            } else {
                string tmp = scratch.file(genomeFileStem(e.name) + ".hgn");
                rewriteHgn(f.f, e.name, e.parent, tmp);
                io = lmap_io_open_file(tmp.c_str());
                _log.push_back("rewrote the genome file of " + e.name);
            }
            e.hgnCrc = put(e.hgnPath, io, &e.hgnBytes);
            lmap_io_close(io);
        }
        if (!e.dnaPath.empty()) {
            auto wasDna = before.find(e.dnaPath);
            auto fromDna = before.find(n.dna.member);
            if (n.dna.store == _in.get() && n.dna.member == e.dnaPath && wasDna != before.end()) {
                e.dnaBytes = wasDna->second->dnaBytes;
                e.dnaCrc = wasDna->second->dnaCrc;
            } else if (w && n.dna.store == _in.get() && fromDna != before.end()) {
                w->link(n.dna.member, e.dnaPath);
                e.dnaBytes = fromDna->second->dnaBytes;
                e.dnaCrc = fromDna->second->dnaCrc;
            } else {
                lmap_io *io = n.dna.store->openMember(n.dna.member);
                e.dnaCrc = put(e.dnaPath, io, &e.dnaBytes);
                lmap_io_close(io);
            }
        }
    }
    if (w) {
        for (const string &member : w->members())
            if (!finalMembers.count(member)) w->remove(member);
        w->putBytes(MANIFEST_NAME, m.text());
        w->commit();
        return;
    }
    /* a directory: new files into place, then the manifest, then what is no longer named */
    for (auto &r : renames)
        if (::rename(r.first.c_str(), r.second.c_str()) != 0)
            throw runtime_error("cannot rename " + r.first + " to " + r.second + ": " + strerror(errno));
    m.write(_path + "/" + MANIFEST_NAME);
    for (const string &member : _in->members())
        if (!finalMembers.count(member)) unlink((_path + "/" + member).c_str());
}
}
