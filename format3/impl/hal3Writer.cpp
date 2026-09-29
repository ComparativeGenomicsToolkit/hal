/*
 * HAL format 3 API: writing alignments.  See hal3/hal3Writer.h.
 */
#include "hal3/hal3Writer.h"

#include "hal3Likelihood.h"
#include "hal3Manifest.h"
#include "hal3Store.h"
#include "liftmap.h"
#include "twobit64.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fcntl.h>
#include <functional>
#include <set>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

namespace hal3 {
namespace api {

/* ---------------------------------------------------------------- one genome */

struct GenomeWriter::Impl {
    string dir, name, parent, stem;
    vector<SequenceDef> seqs, pseqs;
    GenomeWriterOptions opts;
    lmap_builder *b = nullptr;
    tb_writer *dna = nullptr;
    uint32_t dnaSeq = 0;          /* the next sequence to get DNA */
    int64_t dnaIn = -1;           /* bases of the current one; -1 between sequences */
    bool finished = false;

    void check(int rc) {
        if (rc < 0) throw Error(ErrorKind::Io, stem + ".hgn: " + lmap_builder_error(b));
    }
    void abort() {
        if (b) lmap_builder_abort(b);
        b = nullptr;
        if (dna) tb_writer_close(dna, nullptr, 0);   /* incomplete: removes the partial file */
        dna = nullptr;
    }
};

GenomeWriter::GenomeWriter(const string &dir, const string &name, const string &parent, const vector<SequenceDef> &sequences,
                           const vector<SequenceDef> &parentSequences, const GenomeWriterOptions &opts)
    : _impl(new Impl()) {
    Impl &m = *_impl;
    m.dir = dir;
    m.name = name;
    m.parent = parent;
    m.seqs = sequences;
    m.pseqs = parentSequences;
    m.opts = opts;
    m.stem = dir + "/" + genomeFileStem(name);
    if (name.empty()) throw Error(ErrorKind::Argument, "a genome needs a name");
    if (parent.empty() && !parentSequences.empty())
        throw Error(ErrorKind::Argument, "genome " + name + " has no parent but was given parent sequences");
    for (const vector<SequenceDef> *v : {&m.seqs, &m.pseqs}) {
        set<string> seen;
        for (const SequenceDef &s : *v) {
            if (s.name.empty() || s.length < 0)
                throw Error(ErrorKind::Argument, "genome " + name + ": a sequence with no name or a negative length");
            if (!seen.insert(s.name).second)
                throw Error(ErrorKind::Argument, "genome " + name + ": sequence " + s.name + " is named twice" +
                                                     (v == &m.pseqs ? " among its parent's" : ""));
        }
    }
    string hgn = m.stem + ".hgn";
    m.b = lmap_builder_open(hgn.c_str(), HGN_PROFILE, opts.memBytes, opts.tmpDir.empty() ? nullptr : opts.tmpDir.c_str());
    if (!m.b) throw Error(ErrorKind::Io, "cannot write " + hgn + ": " + strerror(errno));
    m.check(lmap_builder_set_params(m.b, opts.chunkRuns, opts.chunkParentSpan, LMAP_CODEC_DEFLATE));
    m.check(lmap_builder_set_order(m.b, LMAP_ORDER_A));
    if (opts.presorted) m.check(lmap_builder_set_presorted(m.b, 1));
    m.check(lmap_builder_set_meta(m.b, "axis.a", ("genome " + name).c_str()));
    m.check(lmap_builder_set_meta(m.b, "hal.genome", name.c_str()));
    if (!parent.empty()) {
        m.check(lmap_builder_set_meta(m.b, "axis.b", ("genome " + parent).c_str()));
        m.check(lmap_builder_set_meta(m.b, "hal.parent", parent.c_str()));
    }
    m.check(lmap_builder_set_meta(m.b, "source", opts.source.c_str()));
    for (const SequenceDef &s : m.seqs) m.check(lmap_builder_add_member(m.b, 0, s.name.c_str(), (uint64_t)s.length));
    for (const SequenceDef &s : m.pseqs) m.check(lmap_builder_add_member(m.b, 1, s.name.c_str(), (uint64_t)s.length));
}

GenomeWriter::~GenomeWriter() {
    if (!_impl->finished) _impl->abort();
}

void GenomeWriter::setMeta(const string &key, const string &value) {
    Impl &m = *_impl;
    if (!m.b) throw Error(ErrorKind::Argument, "genome " + m.name + " is finished");
    m.check(lmap_builder_set_meta(m.b, ("hal.meta." + key).c_str(), value.c_str()));
}

void GenomeWriter::addRun(uint32_t childSeq, int64_t childStart, uint32_t parentSeq, int64_t parentStart, int64_t length,
                          bool reversed) {
    Impl &m = *_impl;
    if (!m.b) throw Error(ErrorKind::Argument, "genome " + m.name + " is finished");
    if (childSeq >= m.seqs.size() || parentSeq >= m.pseqs.size() || length <= 0 || childStart < 0 || parentStart < 0 ||
        childStart + length > m.seqs[childSeq].length || parentStart + length > m.pseqs[parentSeq].length)
        throw Error(ErrorKind::Argument, "genome " + m.name + ": a run outside its sequences (child " + to_string(childSeq) + ":" +
                                             to_string(childStart) + ", parent " + to_string(parentSeq) + ":" +
                                             to_string(parentStart) + ", length " + to_string(length) + ")");
    m.check(lmap_builder_add_run(m.b, childSeq, parentSeq, childStart, parentStart, length, reversed ? 1 : 0));
}

void GenomeWriter::beginDna() {
    Impl &m = *_impl;
    if (!m.b) throw Error(ErrorKind::Argument, "genome " + m.name + " is finished");
    if (m.dnaIn >= 0) throw Error(ErrorKind::Argument, "genome " + m.name + ": beginDna twice");
    if (m.dnaSeq >= m.seqs.size()) throw Error(ErrorKind::Argument, "genome " + m.name + ": DNA for more sequences than it has");
    if (!m.dna) {
        vector<const char *> names;
        for (const SequenceDef &s : m.seqs) {
            if (s.name.size() > 255)
                throw Error(ErrorKind::Argument, "genome " + m.name + ": sequence name longer than a .2bit allows: " + s.name);
            names.push_back(s.name.c_str());
        }
        string tb = m.stem + ".2bit";
        m.dna = tb_writer_open(tb.c_str(), names.data(), (uint32_t)names.size(), 0);
        if (!m.dna) throw Error(ErrorKind::Io, "cannot write " + tb + ": " + strerror(errno));
    }
    if (tb_writer_begin(m.dna) != 0) throw Error(ErrorKind::Io, m.stem + ".2bit: " + tb_writer_error(m.dna));
    m.dnaIn = 0;
}

void GenomeWriter::appendDna(const char *bases, size_t n) {
    Impl &m = *_impl;
    if (m.dnaIn < 0) throw Error(ErrorKind::Argument, "genome " + m.name + ": appendDna outside beginDna / endDna");
    if (m.dnaIn + (int64_t)n > m.seqs[m.dnaSeq].length)
        throw Error(ErrorKind::Argument, "genome " + m.name + ": more DNA than sequence " + m.seqs[m.dnaSeq].name + " has bases");
    if (n && tb_writer_append(m.dna, bases, n) != 0) throw Error(ErrorKind::Io, m.stem + ".2bit: " + tb_writer_error(m.dna));
    m.dnaIn += (int64_t)n;
}

void GenomeWriter::endDna() {
    Impl &m = *_impl;
    if (m.dnaIn < 0) throw Error(ErrorKind::Argument, "genome " + m.name + ": endDna without beginDna");
    if (m.dnaIn != m.seqs[m.dnaSeq].length)
        throw Error(ErrorKind::Argument, "genome " + m.name + ": sequence " + m.seqs[m.dnaSeq].name + " got " + to_string(m.dnaIn) +
                                             " bases of DNA, not " + to_string(m.seqs[m.dnaSeq].length));
    if (tb_writer_end(m.dna) != 0) throw Error(ErrorKind::Io, m.stem + ".2bit: " + tb_writer_error(m.dna));
    m.dnaIn = -1;
    m.dnaSeq++;
}

void GenomeWriter::addDna(const char *bases, int64_t length) {
    beginDna();
    appendDna(bases, (size_t)length);
    endDna();
}

GenomeFiles GenomeWriter::finish() {
    Impl &m = *_impl;
    if (!m.b) throw Error(ErrorKind::Argument, "genome " + m.name + " is finished");
    if (m.dna && (m.dnaIn >= 0 || m.dnaSeq != m.seqs.size())) {
        m.abort();
        throw Error(ErrorKind::Format, "genome " + m.name + ": DNA for " + to_string(m.dnaSeq) + " of its " + to_string(m.seqs.size()) +
                                           " sequences");
    }
    GenomeFiles out;
    out.name = m.name;
    out.parent = m.parent;
    out.numSequences = m.seqs.size();
    for (const SequenceDef &s : m.seqs) out.length += s.length;
    lmap_build_stats st;
    lmap_builder *b = m.b;
    m.b = nullptr;
    if (lmap_builder_close(b, &st) != 0) {
        m.abort();
        throw Error(ErrorKind::Format, "genome " + m.name + ": " + st.error);
    }
    out.runs = st.runs;
    if (m.parent.empty() && st.runs) {
        m.abort();
        unlink((m.stem + ".hgn").c_str());
        throw Error(ErrorKind::Format, "genome " + m.name + " has no parent but has runs");
    }
    out.hgnFile = m.stem + ".hgn";
    if (m.dna) {
        char err[256];
        tb_writer *w = m.dna;
        m.dna = nullptr;
        if (tb_writer_close(w, err, sizeof err) != 0) {
            unlink(out.hgnFile.c_str());
            throw Error(ErrorKind::Io, m.stem + ".2bit: " + err);
        }
        out.dnaFile = m.stem + ".2bit";
    }
    m.finished = true;
    struct stat sb;
    if (m.opts.checksums) {
        out.hgnCrc = fileCrc32(out.hgnFile, &out.hgnBytes);
        if (!out.dnaFile.empty()) out.dnaCrc = fileCrc32(out.dnaFile, &out.dnaBytes);
    } else {
        if (stat(out.hgnFile.c_str(), &sb) == 0) out.hgnBytes = (uint64_t)sb.st_size;
        if (!out.dnaFile.empty() && stat(out.dnaFile.c_str(), &sb) == 0) out.dnaBytes = (uint64_t)sb.st_size;
    }
    return out;
}

/* ---------------------------------------------------------------- manifests */

namespace {

/* The manifest of genomes whose parents are among them: one root, the tree from the
 * parent links (children in the order given), the roster in the tree's preorder. */
Manifest assemble(const vector<GenomeEntry> &entries, const map<string, string> &meta) {
    map<string, size_t> at;
    for (size_t i = 0; i < entries.size(); i++)
        if (!at.emplace(entries[i].name, i).second) throw Error(ErrorKind::Format, "genome " + entries[i].name + " is given twice");
    map<string, vector<string>> kids;
    vector<string> roots;
    for (const GenomeEntry &e : entries) {
        if (e.parent.empty()) roots.push_back(e.name);
        else if (!at.count(e.parent)) throw Error(ErrorKind::Format, "genome " + e.name + "'s parent " + e.parent + " is not here");
        else kids[e.parent].push_back(e.name);
    }
    if (roots.size() != 1) {
        string r;
        for (const string &x : roots) r += " " + x;
        throw Error(ErrorKind::Format, roots.empty() ? "no root (a cycle of parents)" : "more than one root:" + r);
    }
    Manifest m;
    m.meta = meta;
    set<string> placed;
    function<string(const string &)> walk = [&](const string &g) -> string {
        if (!placed.insert(g).second) throw Error(ErrorKind::Format, "genome " + g + " is in a cycle of parents");
        const GenomeEntry &e = entries[at[g]];
        m.genomes.push_back(e);
        string s;
        auto k = kids.find(g);
        if (k != kids.end()) {
            s = "(";
            for (size_t i = 0; i < k->second.size(); i++) s += (i ? "," : "") + walk(k->second[i]);
            s += ")";
        }
        s += quoteNewick(g);
        if (!e.parent.empty() && e.branchLength >= 0) s += ":" + fmtLength(e.branchLength);
        return s;
    };
    m.newick = walk(roots[0]) + ";";
    if (placed.size() != entries.size()) throw Error(ErrorKind::Format, "some genomes are in a cycle of parents");
    return m;
}

bool exists(const string &p) {
    struct stat sb;
    return stat(p.c_str(), &sb) == 0;
}

void mkdirOrThrow(const string &d) {
    if (mkdir(d.c_str(), 0777) != 0) throw Error(ErrorKind::Io, "cannot create " + d + ": " + strerror(errno));
}
}

/* ---------------------------------------------------------------- fragments */

struct Fragment::Impl {
    string dir;
    vector<GenomeEntry> entries;
    map<string, string> meta;
    bool finished = false;
};

Fragment::Fragment(const string &dir) : _impl(new Impl()) {
    _impl->dir = dir;
    if (exists(dir + "/" + MANIFEST_NAME)) throw Error(ErrorKind::Argument, dir + " already holds an alignment");
    if (!exists(dir)) mkdirOrThrow(dir);
    if (!exists(dir + "/genomes")) mkdirOrThrow(dir + "/genomes");
}

Fragment::~Fragment() {}

string Fragment::genomesDir() const { return _impl->dir + "/genomes"; }

void Fragment::add(const GenomeFiles &g, double branchLength) {
    const string gd = genomesDir() + "/";
    auto member = [&](const string &file) {
        if (file.compare(0, gd.size(), gd) != 0)
            throw Error(ErrorKind::Argument, file + " is not in the fragment's genomes directory " + gd);
        return "genomes/" + file.substr(gd.size());
    };
    GenomeEntry e;
    e.name = g.name;
    e.parent = g.parent;
    e.branchLength = g.parent.empty() ? -1 : branchLength;
    e.length = (uint64_t)g.length;
    e.numSequences = g.numSequences;
    e.hgnPath = member(g.hgnFile);
    e.hgnBytes = g.hgnBytes;
    e.hgnCrc = g.hgnCrc;
    if (!g.dnaFile.empty()) {
        e.dnaPath = member(g.dnaFile);
        e.dnaBytes = g.dnaBytes;
        e.dnaCrc = g.dnaCrc;
    }
    _impl->entries.push_back(e);
}

void Fragment::setMeta(const string &key, const string &value) { _impl->meta[key] = value; }

void Fragment::finish() {
    Impl &f = *_impl;
    Manifest m = assemble(f.entries, f.meta);
    for (GenomeEntry &e : m.genomes) {   /* checksums, where the writers did not take them */
        if (!e.hgnCrc) e.hgnCrc = fileCrc32(f.dir + "/" + e.hgnPath, &e.hgnBytes);
        if (!e.dnaPath.empty() && !e.dnaCrc) e.dnaCrc = fileCrc32(f.dir + "/" + e.dnaPath, &e.dnaBytes);
    }
    m.write(f.dir + "/" + MANIFEST_NAME);
    ValidateOptions vo;
    vo.checksums = false;
    vo.deep = false;
    vector<string> problems = validate(f.dir, vo);
    if (!problems.empty()) {
        string s;
        for (size_t i = 0; i < problems.size() && i < 5; i++) s += "; " + problems[i];
        throw Error(ErrorKind::Format, "fragment " + f.dir + ": " + to_string(problems.size()) + " problems" + s);
    }
    f.finished = true;
}

/* ---------------------------------------------------------------- merge */

namespace {

/* One genome as it appears in one fragment. */
struct Seen {
    size_t frag;
    const GenomeEntry *e;
};

/* hard link, else copy (checking the CRC); returns the CRC */
uint32_t placeFile(const Store &from, const string &member, const string &fromDir, const string &path, bool link, uint64_t *bytes) {
    if (link && !fromDir.empty() && ::link((fromDir + "/" + member).c_str(), path.c_str()) == 0)
        return fileCrc32(path, bytes);
    lmap_io *io = from.openMember(member);
    uint32_t c;
    try {
        c = copyToFile(io, path, bytes);
    } catch (...) {
        lmap_io_close(io);
        throw;
    }
    lmap_io_close(io);
    return c;
}
}

void merge(const vector<string> &fragments, const string &out, const map<string, string> &meta, const MergeOptions &opts) {
    if (fragments.empty()) throw Error(ErrorKind::Argument, "nothing to merge");
    if (exists(out)) throw Error(ErrorKind::Argument, out + " already exists");
    vector<unique_ptr<Store>> stores;
    vector<Manifest> manifests;
    for (const string &f : fragments) {
        stores.push_back(Store::open(f));
        manifests.push_back(Manifest::parse(stores.back()->read(MANIFEST_NAME), f));
    }
    /* every genome's appearances, in fragment order; the metadata */
    map<string, vector<Seen>> seen;
    vector<string> order;   /* genomes, first appearance */
    map<string, string> allMeta;
    for (size_t i = 0; i < manifests.size(); i++) {
        for (const GenomeEntry &e : manifests[i].genomes) {
            auto &v = seen[e.name];
            if (v.empty()) order.push_back(e.name);
            v.push_back(Seen{i, &e});
        }
        for (const auto &kv : manifests[i].meta) {
            auto m = allMeta.find(kv.first);
            if (m != allMeta.end() && m->second != kv.second)
                throw Error(ErrorKind::Format, "metadata " + kv.first + " differs between fragments (" + m->second + ", " + kv.second + ")");
            allMeta[kv.first] = kv.second;
        }
    }
    for (const auto &kv : meta) allMeta[kv.first] = kv.second;
    /* each genome: its .hgn where it has a parent, its DNA where it is a root */
    struct Pick {
        GenomeEntry e;
        Seen hgn, dna;
        bool hasDna = false;
    };
    map<string, Pick> picks;
    vector<GenomeEntry> entries;
    for (const string &g : order) {
        const vector<Seen> &v = seen[g];
        const Seen *h = nullptr, *d = nullptr;
        for (const Seen &x : v) {
            if (!x.e->parent.empty()) {
                if (h && (h->e->parent != x.e->parent || h->e->hgnCrc != x.e->hgnCrc))
                    throw Error(ErrorKind::Format, "genome " + g + " has a parent in two fragments, with different branches");
                if (!h) h = &x;
            }
            if (!x.e->dnaPath.empty() && (!d || (x.e->parent.empty() && !d->e->parent.empty()))) d = &x;
            if (x.e->length != v[0].e->length || x.e->numSequences != v[0].e->numSequences)
                throw Error(ErrorKind::Format, "genome " + g + " has different sequences in two fragments");
        }
        if (opts.checkDna)
            for (const Seen &x : v)
                if (d && !x.e->dnaPath.empty() && (x.e->dnaCrc != d->e->dnaCrc || x.e->dnaBytes != d->e->dnaBytes))
                    throw Error(ErrorKind::Format, "genome " + g + "'s DNA differs between fragments " + fragments[d->frag] + " and " +
                                                       fragments[x.frag]);
        if (!h) h = &v[0];
        Pick p;
        p.e = *h->e;
        p.hgn = *h;
        if (d) {
            p.dna = *d;
            p.hasDna = true;
            p.e.dnaBytes = d->e->dnaBytes;
            p.e.dnaCrc = d->e->dnaCrc;
        }
        /* members (e.g. likelihoods): from wherever they are, the same wherever repeated */
        p.e.members.clear();
        map<string, const MemberRef *> roles;
        for (const Seen &x : v)
            for (const MemberRef &mr : x.e->members) {
                auto r = roles.find(mr.role);
                if (r != roles.end() && r->second->crc != mr.crc)
                    throw Error(ErrorKind::Format, "genome " + g + "'s " + mr.role + " member differs between fragments");
                if (r == roles.end()) {
                    roles[mr.role] = &mr;
                    p.e.members.push_back(mr);
                }
            }
        const string stem = "genomes/" + genomeFileStem(g);
        p.e.hgnPath = stem + ".hgn";
        p.e.dnaPath = p.hasDna ? stem + ".2bit" : "";
        for (MemberRef &mr : p.e.members) mr.path = stem + mr.path.substr(mr.path.rfind('/') + 1).substr(genomeFileStem(g).size());
        entries.push_back(p.e);
        picks[g] = p;
    }
    /* siblings in the order of the fragment that holds their branches (where they have a
     * parent), whatever order the fragments come in */
    auto rank = [&](const GenomeEntry &e) {
        const Seen &s = picks[e.name].hgn;
        return make_pair(s.frag, (size_t)(s.e - manifests[s.frag].genomes.data()));
    };
    stable_sort(entries.begin(), entries.end(), [&](const GenomeEntry &x, const GenomeEntry &y) { return rank(x) < rank(y); });
    Manifest m = assemble(entries, allMeta);

    /* the files, each once: into a pack, or a directory */
    auto sourceOf = [&](const Pick &p, int what, const MemberRef *mr, size_t &frag, string &member, uint64_t &bytes, uint32_t &crc) {
        const Seen &s = what == 0 ? p.hgn : what == 1 ? p.dna : p.hgn;
        frag = s.frag;
        if (what == 0) {
            member = s.e->hgnPath;
            bytes = s.e->hgnBytes;
            crc = s.e->hgnCrc;
        } else if (what == 1) {
            member = s.e->dnaPath;
            bytes = s.e->dnaBytes;
            crc = s.e->dnaCrc;
        } else {
            /* a member: from the fragment that has it */
            for (const Seen &x : seen[p.e.name])
                for (const MemberRef &r : x.e->members)
                    if (r.role == mr->role) {
                        frag = x.frag;
                        member = r.path;
                        bytes = r.bytes;
                        crc = r.crc;
                        return;
                    }
        }
    };
    vector<string> made;
    try {
        unique_ptr<PackWriter> pack;
        if (opts.pack) pack = PackWriter::create(out);
        else {
            mkdirOrThrow(out);
            made.push_back(out);
            mkdirOrThrow(out + "/genomes");
            made.push_back(out + "/genomes");
        }
        for (const GenomeEntry &e : m.genomes) {
            const Pick &p = picks[e.name];
            vector<pair<int, const MemberRef *>> files{{0, nullptr}};
            if (p.hasDna) files.push_back({1, nullptr});
            for (const MemberRef &mr : e.members) files.push_back({2, &mr});
            for (auto &f : files) {
                size_t frag = 0;
                string member;
                uint64_t bytes = 0, got = 0;
                uint32_t crc = 0, c;
                sourceOf(p, f.first, f.second, frag, member, bytes, crc);
                const string to = f.first == 0 ? e.hgnPath : f.first == 1 ? e.dnaPath : f.second->path;
                if (pack) {
                    lmap_io *io = stores[frag]->openMember(member);
                    got = (uint64_t)lmap_io_size(io);
                    try {
                        c = pack->putIo(to, io);
                    } catch (...) {
                        lmap_io_close(io);
                        throw;
                    }
                    lmap_io_close(io);
                } else {
                    const string path = out + "/" + to;
                    made.push_back(path);
                    c = placeFile(*stores[frag], member, stores[frag]->isPack() ? "" : fragments[frag], path, opts.link, &got);
                }
                if (c != crc || got != bytes)
                    throw Error(ErrorKind::Checksum, fragments[frag] + ": " + member + " differs from its manifest's size and CRC");
            }
        }
        if (pack) {
            pack->putBytes(MANIFEST_NAME, m.text());
            pack->commit();
        } else {
            m.write(out + "/" + MANIFEST_NAME);
            made.push_back(out + "/" + MANIFEST_NAME);
        }
        ValidateOptions vo;
        vo.checksums = opts.checksums;
        vo.deep = false;
        vector<string> problems = validate(out, vo);
        if (!problems.empty()) {
            string s;
            for (size_t i = 0; i < problems.size() && i < 5; i++) s += "; " + problems[i];
            throw Error(ErrorKind::Format, "the merged alignment has " + to_string(problems.size()) + " problems" + s);
        }
    } catch (...) {
        if (opts.pack) {
            unlink(out.c_str());
            unlink((out + ".partial").c_str());
        } else {
            for (auto i = made.rbegin(); i != made.rend(); ++i)
                if (unlink(i->c_str()) != 0) rmdir(i->c_str());
        }
        throw;
    }
}

/* ---------------------------------------------------------------- validate */

vector<string> validate(const string &path, const ValidateOptions &opts) {
    vector<string> problems;
    auto fail = [&](const string &what) { problems.push_back(what); };
    Manifest m;
    unique_ptr<Store> store;
    try {
        store = Store::open(path);
        m = Manifest::parse(store->read(MANIFEST_NAME), path);
    } catch (exception &e) {
        fail(e.what());
        return problems;
    }

    /* ---- the tree against the roster */
    try {
        vector<NewickNode> nodes = parseNewick(m.newick);
        set<string> seen;
        for (const NewickNode &n : nodes) {
            if (n.name.empty()) {
                fail("the tree has an unnamed node");
                continue;
            }
            if (!seen.insert(n.name).second) fail("the tree names " + n.name + " twice");
            const GenomeEntry *g = m.find(n.name);
            if (!g) {
                fail("genome " + n.name + " is in the tree but not the roster");
                continue;
            }
            string parent = n.parent < 0 ? "" : nodes[n.parent].name;
            if (g->parent != parent) fail(n.name + ": parent is " + parent + " in the tree, " + g->parent + " in the roster");
            if (n.parent >= 0 && n.branchLength >= 0 && fabs(n.branchLength - g->branchLength) > 1e-9 * max(1.0, fabs(n.branchLength)))
                fail(n.name + ": branch length differs between tree and roster");
        }
        for (const GenomeEntry &g : m.genomes)
            if (!seen.count(g.name)) fail("genome " + g.name + " is in the roster but not the tree");
        if (!m.genomes.empty() && !m.genomes[0].parent.empty()) fail("the first genome listed is not the root");
    } catch (exception &e) {
        fail(e.what());
    }

    /* ---- each genome */
    map<string, vector<pair<string, uint64_t>>> tables;   /* genome -> its axis a */
    set<string> named;
    for (const GenomeEntry &g : m.genomes) {
        string gname = "genome " + g.name;
        named.insert(g.hgnPath);
        if (!g.dnaPath.empty()) named.insert(g.dnaPath);
        for (const MemberRef &mr : g.members) {
            if (!named.insert(mr.path).second) fail(gname + ": its " + mr.role + " member " + mr.path + " is named twice");
            if (!store->has(mr.path)) {
                fail(gname + ": its " + mr.role + " member " + mr.path + " is missing");
                continue;
            }
            if (opts.checksums) {
                try {
                    if ((uint64_t)store->size(mr.path) != mr.bytes || store->crc(mr.path) != mr.crc)
                        fail(gname + ": " + mr.path + " size or CRC differs from the manifest");
                } catch (exception &e) {
                    fail(gname + ": " + e.what());
                }
            }
        }
        if (opts.checksums) {
            try {
                if ((uint64_t)store->size(g.hgnPath) != g.hgnBytes || store->crc(g.hgnPath) != g.hgnCrc)
                    fail(gname + ": " + g.hgnPath + " size or CRC differs from the manifest");
                if (!g.dnaPath.empty() && ((uint64_t)store->size(g.dnaPath) != g.dnaBytes || store->crc(g.dnaPath) != g.dnaCrc))
                    fail(gname + ": " + g.dnaPath + " size or CRC differs from the manifest");
            } catch (exception &e) {
                fail(gname + ": " + e.what());
            }
        }
        if (!store->has(g.hgnPath)) {
            fail(gname + ": " + g.hgnPath + " is missing");
            continue;
        }
        unique_ptr<lmap_file, void (*)(lmap_file *)> f(lmap_open_io(store->openMember(g.hgnPath), 1), lmap_close);
        if (!f) {
            fail(gname + ": cannot open " + g.hgnPath);
            continue;
        }
        lmap_file *lf = f.get();
        if (strcmp(lmap_profile(lf), HGN_PROFILE) != 0) fail(gname + ": profile is " + lmap_profile(lf));
        const char *mg = lmap_meta_get(lf, "hal.genome"), *mp = lmap_meta_get(lf, "hal.parent");
        if (!mg || g.name != mg) fail(gname + ": the .hgn says it is " + (mg ? mg : "(none)"));
        if (g.parent != (mp ? mp : "")) fail(gname + ": the .hgn names parent " + (mp ? mp : "(none)"));
        if (opts.deep && lmap_verify(lf) != 0) fail(gname + ": liftmap integrity check failed");
        if (lmap_a_overlap(lf)) fail(gname + ": axis a overlaps -- some base has more than one parent position");
        lmap_file_stats st;
        lmap_get_stats(lf, &st);
        if (g.parent.empty() && (st.runs || lmap_n_members(lf, 1))) fail(gname + ": the root has runs or parent sequences");
        vector<pair<string, uint64_t>> &t = tables[g.name];
        uint64_t total = 0;
        for (uint32_t i = 0; i < lmap_n_members(lf, 0); i++) {
            const lmap_member *mm = lmap_member_at(lf, 0, i);
            t.push_back({mm->name, mm->length});
            total += mm->length;
        }
        if (t.size() != g.numSequences || total != g.length)
            fail(gname + ": " + to_string(t.size()) + " sequences, " + to_string(total) + " bp; the manifest says " +
                 to_string(g.numSequences) + ", " + to_string(g.length));
        if (!g.dnaPath.empty()) {
            unique_ptr<tb_file, void (*)(tb_file *)> d(nullptr, tb_close);
            if (!store->has(g.dnaPath)) fail(gname + ": " + g.dnaPath + " is missing");
            else {
                lmap_io *io = store->openMember(g.dnaPath);
                d.reset(tb_open_io(io, 1));
                if (!d) {
                    lmap_io_close(io);
                    fail(gname + ": cannot open " + g.dnaPath);
                }
            }
            if (d) {
                if (tb_n_seqs(d.get()) != t.size()) fail(gname + ": the .2bit has " + to_string(tb_n_seqs(d.get())) + " sequences");
                else
                    for (uint32_t i = 0; i < t.size(); i++) {
                        char err[256];
                        if (t[i].first != tb_name(d.get(), i) || tb_length(d.get(), i) != (int64_t)t[i].second) {
                            fail(gname + ": .2bit sequence " + to_string(i) + " is " + tb_name(d.get(), i) + ", not " + t[i].first);
                            break;
                        }
                        if (opts.deep && tb_verify_seq(d.get(), i, err, sizeof err) != 0) fail(gname + ": .2bit " + t[i].first + ": " + err);
                    }
            }
        }
        /* likelihoods: for this genome's sequences, every chunk sound */
        if (const MemberRef *mr = g.member("lik")) {
            if (store->has(mr->path)) {
                try {
                    unique_ptr<LikFile> lk = LikFile::open(store->openMember(mr->path), true, mr->path);
                    bool same = lk->numSequences() == t.size();
                    for (uint32_t i = 0; same && i < t.size(); i++)
                        same = lk->sequenceName(i) == t[i].first && lk->sequenceLength(i) == (int64_t)t[i].second;
                    if (!same) fail(gname + ": " + mr->path + " is not for its sequences");
                    else if (opts.deep)
                        for (const string &p : lk->verify()) fail(gname + ": " + p);
                    if (g.dnaPath.empty()) fail(gname + ": has likelihoods but no DNA (they are relative to its bases)");
                } catch (exception &e) {
                    fail(gname + ": " + e.what());
                }
            }
        }
        /* axis b against the parent's axis a (the parent comes first in the roster) */
        if (!g.parent.empty()) {
            auto pt = tables.find(g.parent);
            if (pt == tables.end()) fail(gname + ": listed before its parent " + g.parent);
            else {
                uint32_t nb = lmap_n_members(lf, 1);
                if (nb != pt->second.size())
                    fail(gname + ": axis b has " + to_string(nb) + " sequences, the parent " + to_string(pt->second.size()));
                else
                    for (uint32_t i = 0; i < nb; i++) {
                        const lmap_member *mm = lmap_member_at(lf, 1, i);
                        if (pt->second[i].first != mm->name || pt->second[i].second != mm->length) {
                            fail(gname + ": parent sequence " + mm->name + " differs from " + g.parent + "'s own table");
                            break;
                        }
                    }
            }
        }
    }

    /* ---- members the manifest does not name; in a pack, members against its index */
    for (const string &n : store->members())
        if (n != MANIFEST_NAME && !named.count(n)) fail(n + " is not in the manifest");
    if (store->isPack() && opts.checksums)
        for (const auto &e : packInfo(*store).entries)
            if (store->crc(e.first) != e.second.crc) fail(e.first + ": CRC differs from the pack's index");
    return problems;
}
}
}
