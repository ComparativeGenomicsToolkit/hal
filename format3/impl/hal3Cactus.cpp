/*
 * HAL format 3 API: a Cactus subproblem (.c2h form) to a fragment.  See
 * hal3/hal3Cactus.h.
 */
#include "hal3/hal3Cactus.h"

#include "hal3Manifest.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

using namespace std;

namespace hal3 {
namespace api {

namespace {

/* the fields of one tab-separated line [p, e) */
vector<const char *> fields(const char *p, const char *e, vector<const char *> &ends) {
    vector<const char *> f;
    ends.clear();
    const char *s = p;
    for (const char *q = p;; q++) {
        if (q == e || *q == '\t') {
            f.push_back(s);
            ends.push_back(q);
            if (q == e) break;
            s = q + 1;
        }
    }
    return f;
}

int64_t number(const char *s, const char *e, const string &what) {
    if (s == e) throw Error(ErrorKind::Format, "an empty number in " + what);
    int64_t v = 0;
    bool neg = *s == '-';
    for (const char *q = neg ? s + 1 : s; q < e; q++) {
        if (*q < '0' || *q > '9') throw Error(ErrorKind::Format, "not a number (" + string(s, e) + ") in " + what);
        v = v * 10 + (*q - '0');
    }
    return neg ? -v : v;
}
}

struct C2hImporter::Impl {
    unique_ptr<Fragment> fragment;
    GenomeWriterOptions opts;
    string root;
    vector<string> children;             /* in the tree's order */
    map<string, double> branch;          /* child -> branch length */
    struct Genome {
        vector<SequenceDef> seqs;
        map<string, uint32_t> index;
        uint32_t next = 0;               /* the next sequence the second pass expects */
        unique_ptr<GenomeWriter> writer;
    };
    map<string, Genome> genomes;         /* kept genomes: the root and its children */
    struct Bottom {
        uint32_t seq;
        int64_t start, length;
    };
    unordered_map<int64_t, Bottom> bottoms;   /* the ancestor's segments, by name */
    bool finished = false;

    bool kept(const string &g) const { return genomes.count(g) != 0; }
};

C2hImporter::C2hImporter(const string &fragmentDir, const string &newick, const GenomeWriterOptions &opts) : _impl(new Impl()) {
    Impl &m = *_impl;
    m.opts = opts;
    vector<NewickNode> nodes = parseNewick(newick);
    if (nodes.empty() || nodes[0].name.empty()) throw Error(ErrorKind::Format, "the subproblem's tree has no named root: " + newick);
    m.root = nodes[0].name;
    m.genomes[m.root];
    for (const NewickNode &n : nodes) {
        if (n.parent != 0) {
            if (n.parent > 0) throw Error(ErrorKind::Format, "the subproblem's tree is deeper than an ancestor and its children: " + newick);
            continue;
        }
        if (n.name.empty()) throw Error(ErrorKind::Format, "the subproblem's tree has an unnamed child: " + newick);
        m.children.push_back(n.name);
        m.branch[n.name] = n.branchLength;
        m.genomes[n.name];
    }
    m.fragment.reset(new Fragment(fragmentDir));
}

C2hImporter::~C2hImporter() {}

void C2hImporter::declare(const string &genome, const string &name, bool ancestor, int64_t length) {
    Impl &m = *_impl;
    if (!m.kept(genome)) return;   /* an outgroup */
    if (ancestor != (genome == m.root))
        throw Error(ErrorKind::Format, "sequence " + name + " of " + genome + (ancestor ? " is" : " is not") +
                                           " marked as the ancestor's, but the subproblem's ancestor is " + m.root);
    Impl::Genome &g = m.genomes[genome];
    if (!g.index.emplace(name, (uint32_t)g.seqs.size()).second)
        throw Error(ErrorKind::Format, "genome " + genome + " has two sequences named " + name);
    g.seqs.push_back(SequenceDef{name, length});
}

void C2hImporter::sequence(const string &genome, const string &name, const char *segments, size_t segmentsLength, const char *dna,
                           int64_t dnaLength) {
    Impl &m = *_impl;
    if (!m.kept(genome)) return;
    Impl::Genome &g = m.genomes[genome];
    if (g.next >= g.seqs.size() || g.seqs[g.next].name != name)
        throw Error(ErrorKind::Format, "sequence " + name + " of " + genome + " comes out of the first pass's order");
    const uint32_t si = g.next++;
    const bool isRoot = genome == m.root;
    if (!g.writer) {
        if (isRoot) {
            g.writer.reset(new GenomeWriter(m.fragment->genomesDir(), genome, "", g.seqs, {}, m.opts));
        } else {
            const Impl::Genome &r = m.genomes[m.root];
            if (r.writer == nullptr && !r.seqs.empty() && r.next < r.seqs.size())
                throw Error(ErrorKind::Format, "genome " + genome + "'s segments come before all of the ancestor's");
            g.writer.reset(new GenomeWriter(m.fragment->genomesDir(), genome, m.root, g.seqs, r.seqs, m.opts));
        }
    }
    const int64_t len = g.seqs[si].length;
    const string what = "sequence " + name + " of " + genome;
    vector<const char *> f, ends;
    int64_t covered = 0;
    const char *p = segments, *end = segments + segmentsLength;
    while (p < end) {
        const char *e = (const char *)memchr(p, '\n', (size_t)(end - p));
        if (!e) e = end;
        if (e > p) {
            if (*p != 'a') throw Error(ErrorKind::Format, "a line that is not a segment in " + what + ": " + string(p, e));
            f = fields(p, e, ends);
            if (isRoot) {   /* a <name> <start> <length> */
                if (f.size() != 4) throw Error(ErrorKind::Format, "an ancestor segment without 4 fields in " + what);
                int64_t nm = number(f[1], ends[1], what), s = number(f[2], ends[2], what), l = number(f[3], ends[3], what);
                if (s < 0 || l <= 0 || s + l > len) throw Error(ErrorKind::Format, "a segment outside " + what);
                if (!m.bottoms.emplace(nm, Impl::Bottom{si, s, l}).second)
                    throw Error(ErrorKind::Format, "two ancestor segments named " + to_string(nm));
                covered += l;
            } else {        /* a <start> <length> [<ancestor segment> <reversed>] */
                if (f.size() != 3 && f.size() != 5) throw Error(ErrorKind::Format, "a segment without 3 or 5 fields in " + what);
                int64_t s = number(f[1], ends[1], what), l = number(f[2], ends[2], what);
                if (s < 0 || l <= 0 || s + l > len) throw Error(ErrorKind::Format, "a segment outside " + what);
                covered += l;
                if (f.size() == 5) {
                    int64_t nm = number(f[3], ends[3], what), rv = number(f[4], ends[4], what);
                    auto b = m.bottoms.find(nm);
                    if (b == m.bottoms.end())
                        throw Error(ErrorKind::Format, what + " is aligned to ancestor segment " + to_string(nm) + ", which it does not have");
                    if (b->second.length != l)
                        throw Error(ErrorKind::Format, what + ": a segment of " + to_string(l) + " bases aligned to an ancestor segment of " +
                                                           to_string(b->second.length));
                    /* the field is the ancestor segment's strand: 1 the same as the child's */
                    g.writer->addRun(si, s, b->second.seq, b->second.start, l, rv == 0);
                }
            }
        }
        p = e + 1;
    }
    if (covered != len)
        throw Error(ErrorKind::Format, what + ": its segments cover " + to_string(covered) + " of its " + to_string(len) + " bases");
    g.writer->addDna(dna, dnaLength);
}

void C2hImporter::finish() {
    Impl &m = *_impl;
    if (m.finished) throw Error(ErrorKind::Argument, "the subproblem is finished");
    vector<string> order{m.root};
    order.insert(order.end(), m.children.begin(), m.children.end());
    for (const string &name : order) {
        Impl::Genome &g = m.genomes[name];
        if (g.seqs.empty()) throw Error(ErrorKind::Format, "genome " + name + " of the subproblem's tree has no sequences in it");
        if (g.next != g.seqs.size())
            throw Error(ErrorKind::Format, "genome " + name + ": " + to_string(g.next) + " of its " + to_string(g.seqs.size()) +
                                               " sequences given");
        m.fragment->add(g.writer->finish(), name == m.root ? -1 : m.branch[name]);
    }
    m.fragment->finish();
    m.finished = true;
}

/* ---------------------------------------------------------------- from files */

namespace {

/* A FASTA file's records, by name, read on demand. */
class Fasta {
  public:
    explicit Fasta(const string &path) : _path(path), _in(path, ios::binary) {
        if (!_in) throw Error(ErrorKind::Io, "cannot read " + path);
        string line;
        streamoff at = 0;
        while (getline(_in, line)) {
            streamoff next = at + (streamoff)line.size() + 1;
            if (!line.empty() && line[0] == '>') {
                size_t e = line.find_first_of(" \t\r", 1);
                string name = line.substr(1, e == string::npos ? string::npos : e - 1);
                if (!_at.emplace(name, next).second) throw Error(ErrorKind::Format, path + " has two records named " + name);
            }
            at = next;
        }
        _in.clear();
    }
    /* the record's bases */
    void get(const string &name, string &out) {
        auto a = _at.find(name);
        if (a == _at.end()) throw Error(ErrorKind::Format, _path + " has no record " + name);
        out.clear();
        _in.clear();
        _in.seekg(a->second);
        string line;
        while (getline(_in, line) && (line.empty() || line[0] != '>')) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            out += line;
        }
    }

  private:
    string _path;
    ifstream _in;
    map<string, streamoff> _at;
};

/* each sequence of a .c2h: its "s" line's fields and its segment lines */
void scanC2h(const string &path, const function<void(const string &, const string &, bool, const string &)> &visit) {
    ifstream in(path, ios::binary);
    if (!in) throw Error(ErrorKind::Io, "cannot read " + path);
    string line, genome, name, segs;
    bool ancestor = false, have = false;
    auto unquote = [&](const string &s) { return s.size() >= 2 && s.front() == '\'' && s.back() == '\'' ? s.substr(1, s.size() - 2) : s; };
    while (getline(in, line)) {
        if (!line.empty() && line[0] == 's') {
            if (have) visit(genome, name, ancestor, segs);
            vector<string> f;
            stringstream ss(line);
            string t;
            while (getline(ss, t, '\t')) f.push_back(t);
            if (f.size() != 4) throw Error(ErrorKind::Format, path + ": a sequence line without 4 fields: " + line);
            genome = unquote(f[1]);
            name = unquote(f[2]);
            ancestor = f[3] == "1";
            segs.clear();
            have = true;
        } else if (!line.empty()) {
            if (!have) throw Error(ErrorKind::Format, path + ": a segment before any sequence");
            segs += line;
            segs += '\n';
        }
    }
    if (have) visit(genome, name, ancestor, segs);
}

/* the total length of a sequence's segment lines (the second numeric field of a child
 * segment, the third of an ancestor's) */
int64_t segmentsLength(const string &segs, bool ancestor, const string &what) {
    int64_t total = 0;
    vector<const char *> f, ends;
    const char *p = segs.data(), *end = p + segs.size();
    while (p < end) {
        const char *e = (const char *)memchr(p, '\n', (size_t)(end - p));
        if (!e) e = end;
        if (e > p) {
            f = fields(p, e, ends);
            size_t k = ancestor ? 3 : 2;
            if (f.size() <= k) throw Error(ErrorKind::Format, "a short segment line in " + what);
            total += number(f[k], ends[k], what);
        }
        p = e + 1;
    }
    return total;
}
}

void importC2h(const string &c2hPath, const string &fastaPath, const string &newick, const string &fragmentDir,
               const GenomeWriterOptions &opts) {
    string tree = newick;
    {
        ifstream t(newick);
        if (t && newick.find('(') == string::npos) {
            stringstream ss;
            ss << t.rdbuf();
            tree = ss.str();
            while (!tree.empty() && isspace((unsigned char)tree.back())) tree.pop_back();
        }
    }
    C2hImporter im(fragmentDir, tree, opts);
    scanC2h(c2hPath, [&](const string &g, const string &n, bool anc, const string &segs) {
        im.declare(g, n, anc, segmentsLength(segs, anc, "sequence " + n + " of " + g));
    });
    Fasta fa(fastaPath);
    string dna;
    set<string> kept;
    for (const NewickNode &n : parseNewick(tree))
        if (!n.name.empty()) kept.insert(n.name);
    scanC2h(c2hPath, [&](const string &g, const string &n, bool, const string &segs) {
        if (!kept.count(g)) return;   /* an outgroup: its DNA is not needed */
        fa.get(n, dna);
        im.sequence(g, n, segs.data(), segs.size(), dna.data(), (int64_t)dna.size());
    });
    im.finish();
}
}
}
