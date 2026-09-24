/*
 * HAL format 3: MAF export, row by row.  See hal3MafExport.h.
 */
#include "hal3MafExport.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <zlib.h>

using namespace std;

namespace hal3 {

/* ---------------------------------------------------------------- output */

namespace {

class PlainSink : public MafSink {
  public:
    explicit PlainSink(const string &path) : _own(path != "stdout") {
        _f = _own ? fopen(path.c_str(), "w") : stdout;
        if (!_f) throw runtime_error("cannot write " + path);
    }
    ~PlainSink() override {
        if (_f) close();
    }
    void write(const string &text) override {
        if (fwrite(text.data(), 1, text.size(), _f) != text.size()) throw runtime_error("write failed");
    }
    void close() override {
        if (!_f) return;
        bool bad = _own ? fclose(_f) != 0 : fflush(_f) != 0;
        _f = nullptr;
        if (bad) throw runtime_error("write failed");
    }

  private:
    FILE *_f;
    bool _own;
};

/* BGZF: gzip members of at most 64 KiB, each with a BC extra field giving its size,
 * ending with the standard empty block. */
class BgzfSink : public MafSink {
  public:
    explicit BgzfSink(const string &path) : _plain(path) {}
    ~BgzfSink() override {
        try {
            close();
        } catch (...) {
        }
    }
    void write(const string &text) override {
        _buf += text;
        while (_buf.size() >= BLOCK) {
            block(_buf.data(), BLOCK);
            _buf.erase(0, BLOCK);
        }
    }
    void close() override {
        if (_closed) return;
        _closed = true;
        if (!_buf.empty()) block(_buf.data(), _buf.size());
        _buf.clear();
        static const unsigned char eof[28] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43,
                                              2,    0,    0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        _plain.write(string((const char *)eof, sizeof eof));
        _plain.close();
    }

  private:
    static const size_t BLOCK = 65280;
    PlainSink _plain;
    string _buf;
    bool _closed = false;
    void block(const char *data, size_t n) {
        unsigned char out[65536];
        z_stream z;
        memset(&z, 0, sizeof z);
        if (deflateInit2(&z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
            throw runtime_error("deflateInit2 failed");
        z.next_in = (Bytef *)data;
        z.avail_in = (uInt)n;
        z.next_out = out + 18;
        z.avail_out = sizeof out - 18 - 8;
        int rc = deflate(&z, Z_FINISH);
        size_t clen = z.total_out;
        deflateEnd(&z);
        if (rc != Z_STREAM_END) throw runtime_error("BGZF block did not fit");
        size_t total = 18 + clen + 8;
        unsigned char hdr[18] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0,
                                 (unsigned char)((total - 1) & 0xff), (unsigned char)((total - 1) >> 8)};
        memcpy(out, hdr, 18);
        uLong crc = crc32(crc32(0L, Z_NULL, 0), (const Bytef *)data, (uInt)n);
        unsigned char *t = out + 18 + clen;
        for (int k = 0; k < 4; k++) t[k] = (unsigned char)(crc >> (8 * k));
        for (int k = 0; k < 4; k++) t[4 + k] = (unsigned char)(n >> (8 * k));
        _plain.write(string((const char *)out, total));
    }
};
}

unique_ptr<MafSink> openMafSink(const string &path, bool bgzf) {
    if (bgzf) return unique_ptr<MafSink>(new BgzfSink(path));
    return unique_ptr<MafSink>(new PlainSink(path));
}

/* ---------------------------------------------------------------- scope */

static WalkOptions walkOptions(const MafOptions &o) {
    WalkOptions w;
    w.ancestors = o.ancestors;
    w.dupes = o.dupes;
    w.onlyOrthologs = o.onlyOrthologs;
    w.targets = o.targets;
    w.rootGenome = o.rootGenome;
    return w;
}

MafExporter::MafExporter(const Alignment &alignment, const Genome *ref, const MafOptions &opts, MafSink &out)
    : _aln(alignment), _ref(ref), _opts(opts), _out(out), _walk(alignment, ref, walkOptions(opts)) {}

MafExporter::~MafExporter() {
    for (auto &d : _dna)
        if (d.second) tb_close(d.second);
}

void MafExporter::writeHeader() { _out.write("##maf version=1 scoring=N/A\n# hal " + _aln.manifest().newick + "\n\n"); }

void MafExporter::exportInterval(uint32_t refSeq, int64_t start, int64_t end) {
    for (int64_t s = start; s < end; s += _opts.chunk) chunk(refSeq, s, min(end, s + _opts.chunk));
}

/* ---------------------------------------------------------------- insertions */

bool MafExporter::Key::operator<(const Key &o) const {
    if (node != o.node) return node->name() < o.node->name();
    if (seq != o.seq) return seq < o.seq;
    return pos < o.pos;
}

/* The bases [o1, o2) of g's sequence seq, in oriented coordinates (on the reverse strand
 * when reversed), as candidate inserted bases.  keys gets each base's top ancestral base,
 * in oriented order.  Returns false if any of them is in a column holding a reference
 * base -- such a base belongs in that column, not in an insertion.
 *
 * Default and --noDupes: a base's column is the down-closure of its top base (walking up
 * with the mode's copy rule), so it holds a reference base exactly when walking down from
 * the top towards the reference reaches one.  --onlyOrthologs: when liftover from the
 * base to the reference (through their common ancestor) reaches one. */
bool MafExporter::insertionKeys(const Genome *g, uint32_t seq, bool reversed, int64_t o1, int64_t o2, vector<Key> &keys) {
    int64_t len = o2 - o1, seqLen = g->sequenceLength(seq);
    int64_t lo = reversed ? seqLen - o2 : o1;   /* forward coordinates */
    keys.assign((size_t)len, Key{nullptr, 0, 0});
    if (_opts.onlyOrthologs) {
        vector<LiftPiece> hits;
        liftInterval(g, seq, lo, lo + len, _ref, _opts.dupes, hits, &_walk.cursors());
        if (!hits.empty()) return false;
    }
    /* up to the tops */
    vector<LiftPiece> cur{LiftPiece{lo, seq, lo, len, false}}, next, unmapped;
    vector<pair<const Genome *, vector<LiftPiece>>> tops;
    for (const Genome *a = g; !cur.empty(); a = a->parent()) {
        if (a == _walk.top()) {
            tops.emplace_back(a, cur);
            break;
        }
        liftUpBranch(a, cur, next, &unmapped, _walk.cursors(), !_opts.dupes);
        if (!unmapped.empty()) tops.emplace_back(a, unmapped);
        cur.swap(next);
    }
    for (auto &t : tops) {
        for (const LiftPiece &p : t.second) {
            for (int64_t i = 0; i < p.length; i++) {
                int64_t fwd = p.srcStart + i;                    /* the base in g, forward */
                int64_t at = p.reversed ? p.tgtStart + p.length - 1 - i : p.tgtStart + i;
                int64_t k = reversed ? (lo + len - 1 - fwd) : (fwd - lo);   /* oriented index */
                keys[(size_t)k] = Key{t.first, p.tgtSeq, at};
            }
        }
        if (_opts.onlyOrthologs) continue;
        /* down from the top towards the reference, if the top is on its lineage */
        const vector<const Genome *> &path = _walk.refPath();
        auto on = find(path.begin(), path.end(), t.first);
        if (on == path.end()) continue;
        vector<LiftPiece> down = t.second;
        for (auto c = on; c != path.begin() && !down.empty();) {
            --c;
            liftDownBranch(*c, _opts.dupes, down, next, _walk.cursors());
            down.swap(next);
        }
        if (!down.empty()) return false;
    }
    /* bases with no key cannot happen (every base tops out somewhere) */
    for (const Key &k : keys)
        if (!k.node) throw runtime_error("internal error: an inserted base with no top");
    return true;
}

/* ---------------------------------------------------------------- rows and blocks */

/* A row: pieces of one genome, contiguous in it except for inserted bases between them,
 * in reference order.  Oriented coordinates (MAF's: position on the row's strand)
 * increase along the reference.  Link i joins pieces i and i+1; the bases between them,
 * [oend(i), ostart[i+1]), are inserted: they go in the insertion group at reference
 * boundary at[i] (between reference positions at[i]-1 and at[i]). */
struct MafExporter::Chain {
    const Genome *genome;
    uint32_t seq;
    bool reversed;
    vector<LiftPiece> pieces;   /* in reference order */
    vector<int64_t> ostart;     /* oriented start of each piece */
    vector<int64_t> at;         /* per link: its boundary, or -1 if nothing is inserted */
    vector<vector<Key>> keys;   /* per link: each inserted base's top key */
    vector<vector<int>> cols;   /* per link: each inserted base's column in its group */
    int64_t start() const { return pieces.front().srcStart; }
    int64_t end() const { return pieces.back().srcStart + pieces.back().length; }
    int64_t oend(size_t i) const { return ostart[i] + pieces[i].length; }
};

namespace {
char complement(char c) {
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'a': return 't';
    case 'c': return 'g';
    case 'g': return 'c';
    case 't': return 'a';
    default: return c;
    }
}
}

void MafExporter::fetch(const Genome *g, uint32_t seq, int64_t start, int64_t end, bool reversed, string &out) {
    out.assign((size_t)(end - start), 'N');
    if (_opts.noDna) return;
    auto it = _dna.find(g);
    if (it == _dna.end()) {
        const GenomeEntry *e = _aln.manifest().find(g->name());
        if (!e || e->dnaPath.empty())
            throw runtime_error("genome " + g->name() + " has no DNA in this alignment (--noDna writes N)");
        tb_file *f = tb_open((_aln.dir() + "/" + e->dnaPath).c_str());
        if (!f) throw runtime_error("cannot open " + e->dnaPath);
        it = _dna.emplace(g, f).first;
        vector<int64_t> &idx = _dnaIndex[g];
        for (uint32_t i = 0; i < g->numSequences(); i++)
            idx.push_back(i < tb_n_seqs(f) && !strcmp(tb_name(f, i), g->sequenceName(i)) ? (int64_t)i
                                                                                        : tb_find(f, g->sequenceName(i)));
    }
    int64_t ti = _dnaIndex[g][seq];
    if (ti < 0) throw runtime_error("sequence " + string(g->sequenceName(seq)) + " is not in the .2bit of " + g->name());
    int64_t len = g->sequenceLength(seq);
    int64_t fs = reversed ? len - end : start, fe = reversed ? len - start : end;   /* forward coordinates */
    if (tb_fetch(it->second, (uint32_t)ti, fs, fe, &out[0], 0) != 0)
        throw runtime_error("cannot read DNA of " + g->name() + "." + g->sequenceName(seq));
    if (reversed) {
        reverse(out.begin(), out.end());
        for (char &c : out) c = complement(c);
    }
}

void MafExporter::chunk(uint32_t refSeq, int64_t cs, int64_t ce) {
    map<const Genome *, vector<LiftPiece>> lifted;
    _walk.lift(refSeq, {{cs, ce}}, lifted);
    const int64_t n = max<int64_t>(0, _opts.joinGap);

    /* ---- rows: link each genome's pieces where one follows another in the genome --
     * directly, or across at most n bases that are in no reference column -- and later
     * on the reference.  Links are taken greedily: fewest inserted bases, then fewest
     * skipped reference bases. */
    vector<Chain> chains;
    for (auto &gp : lifted) {
        const Genome *g = gp.first;
        vector<LiftPiece> &ps = gp.second;
        size_t np = ps.size();
        vector<int64_t> os(np);
        for (size_t i = 0; i < np; i++) {
            const LiftPiece &p = ps[i];
            os[i] = p.reversed ? g->sequenceLength(p.tgtSeq) - (p.tgtStart + p.length) : p.tgtStart;
        }
        map<pair<uint32_t, bool>, vector<pair<int64_t, size_t>>> byStart;   /* oriented start -> piece */
        for (size_t i = 0; i < np; i++) byStart[make_pair(ps[i].tgtSeq, ps[i].reversed)].push_back({os[i], i});
        for (auto &b : byStart) sort(b.second.begin(), b.second.end());
        vector<tuple<int64_t, int64_t, int64_t, size_t, size_t>> edges;   /* dTgt, dRef, p.src, p, q */
        for (size_t i = 0; i < np; i++) {
            const vector<pair<int64_t, size_t>> &v = byStart[make_pair(ps[i].tgtSeq, ps[i].reversed)];
            int64_t oe = os[i] + ps[i].length, pend = ps[i].srcStart + ps[i].length;
            for (auto it = lower_bound(v.begin(), v.end(), make_pair(oe, (size_t)0)); it != v.end() && it->first <= oe + n; ++it)
                if (ps[it->second].srcStart >= pend) edges.emplace_back(it->first - oe, ps[it->second].srcStart - pend, ps[i].srcStart, i, it->second);
        }
        sort(edges.begin(), edges.end());
        vector<int64_t> succ(np, -1), pred(np, -1);
        vector<bool> blocked(np, false);   /* p's nearest gap holds a reference-aligned base */
        map<size_t, vector<Key>> linkKeys;  /* by p */
        vector<Key> keys;
        for (auto &e : edges) {
            int64_t dTgt = get<0>(e);
            size_t p = get<3>(e), q = get<4>(e);
            if (succ[p] >= 0 || pred[q] >= 0 || blocked[p]) continue;
            if (dTgt > 0) {
                /* a wider gap for p contains this one, so one refusal settles p */
                if (!insertionKeys(g, ps[p].tgtSeq, ps[p].reversed, os[p] + ps[p].length, os[q], keys)) {
                    blocked[p] = true;
                    _stats.joinsRejected++;
                    continue;
                }
                linkKeys[p] = keys;
            }
            succ[p] = (int64_t)q;
            pred[q] = (int64_t)p;
        }
        for (size_t i = 0; i < np; i++) {
            if (pred[i] >= 0) continue;
            Chain c{g, ps[i].tgtSeq, ps[i].reversed, {}, {}, {}, {}, {}};
            for (int64_t k = (int64_t)i; k >= 0; k = succ[k]) {
                if (!c.pieces.empty()) {
                    auto lk = linkKeys.find((size_t)pred[k]);
                    c.at.push_back(lk == linkKeys.end() ? -1 : 0);
                    c.keys.push_back(lk == linkKeys.end() ? vector<Key>() : lk->second);
                    c.cols.emplace_back();
                }
                c.pieces.push_back(ps[k]);
                c.ostart.push_back(os[k]);
            }
            chains.push_back(move(c));
        }
    }
    sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        if (a.start() != b.start()) return a.start() < b.start();
        int c = a.genome->name().compare(b.genome->name());
        if (c) return c < 0;
        return make_tuple(a.seq, a.ostart[0], a.reversed) < make_tuple(b.seq, b.ostart[0], b.reversed);
    });

    /* ---- insertion groups.  A link that skips no reference base goes at the boundary
     * where its pieces meet; one that also skips reference bases (a deletion as well)
     * goes at whichever end shares more top keys with what is placed there already
     * (ties: the left), with its deletion gaps after it. */
    struct Member {
        size_t chain, link;
    };
    struct Group {
        vector<Member> members;
        set<Key> keys;
        int width = 0;
        int64_t bases = 0;
    };
    map<int64_t, Group> groups;
    auto place = [&](size_t ci, size_t li, int64_t b) {
        Chain &c = chains[ci];
        c.at[li] = b;
        Group &gr = groups[b];
        gr.members.push_back(Member{ci, li});
        gr.keys.insert(c.keys[li].begin(), c.keys[li].end());
        gr.bases += (int64_t)c.keys[li].size();
    };
    for (int pass = 0; pass < 2; pass++) {
        for (size_t ci = 0; ci < chains.size(); ci++) {
            Chain &c = chains[ci];
            for (size_t li = 0; li < c.at.size(); li++) {
                if (c.at[li] < 0) continue;
                int64_t left = c.pieces[li].srcStart + c.pieces[li].length, right = c.pieces[li + 1].srcStart;
                if ((pass == 0) != (left == right)) continue;
                if (left == right) {
                    place(ci, li, left);
                    continue;
                }
                auto shared = [&](int64_t b) {
                    auto g = groups.find(b);
                    if (g == groups.end()) return 0;
                    int k = 0;
                    for (const Key &key : c.keys[li]) k += (int)g->second.keys.count(key);
                    return k;
                };
                place(ci, li, shared(right) > shared(left) ? right : left);
            }
        }
    }
    /* Layout: members in output order; each base joins the first column after the row's
     * previous base that holds its key, or gets a new column right after it. */
    for (auto &gb : groups) {
        Group &gr = gb.second;
        sort(gr.members.begin(), gr.members.end(), [&](const Member &x, const Member &y) {
            const Chain &a = chains[x.chain], &b = chains[y.chain];
            int c = a.genome->name().compare(b.genome->name());
            if (c) return c < 0;
            return make_tuple(a.seq, a.oend(x.link), a.reversed) < make_tuple(b.seq, b.oend(y.link), b.reversed);
        });
        vector<Key> colKey;
        vector<int> order, pos;
        map<Key, vector<int>> withKey;
        for (const Member &m : gr.members) {
            Chain &c = chains[m.chain];
            vector<int> &ids = c.cols[m.link];
            ids.clear();
            int prev = -1;
            for (const Key &k : c.keys[m.link]) {
                int found = -1;
                auto w = withKey.find(k);
                if (w != withKey.end())
                    for (int id : w->second)
                        if (pos[id] > prev && (found < 0 || pos[id] < pos[found])) found = id;
                if (found < 0) {
                    found = (int)colKey.size();
                    colKey.push_back(k);
                    order.insert(order.begin() + (prev + 1), found);
                    pos.resize(colKey.size());
                    for (size_t i = 0; i < order.size(); i++) pos[order[i]] = (int)i;
                    withKey[k].push_back(found);
                }
                ids.push_back(found);
                prev = pos[found];
            }
        }
        for (const Member &m : gr.members)
            for (int &id : chains[m.chain].cols[m.link]) id = pos[id];
        gr.width = (int)order.size();
    }

    /* ---- counts for the block rule: chain starts and ends, covered bases, groups */
    vector<int64_t> starts, ends;
    vector<pair<int64_t, int>> cov;
    for (const Chain &c : chains) {
        starts.push_back(c.start());
        ends.push_back(c.end());
        for (const LiftPiece &p : c.pieces) {
            cov.push_back({p.srcStart, 1});
            cov.push_back({p.srcStart + p.length, -1});
        }
    }
    sort(starts.begin(), starts.end());
    sort(ends.begin(), ends.end());
    sort(cov.begin(), cov.end());
    vector<int64_t> covPos{cs}, covCum{0};   /* aligned reference-column bases before covPos[i] */
    vector<int> covDepth{0};
    {
        int depth = 0;
        for (auto &e : cov) {
            if (e.first != covPos.back()) {
                covCum.push_back(covCum.back() + (int64_t)depth * (e.first - covPos.back()));
                covPos.push_back(e.first);
                covDepth.push_back(depth);
            }
            depth += e.second;
            covDepth.back() = depth;
        }
    }
    auto aligned = [&](int64_t x) {   /* aligned bases in reference columns [cs, x) */
        size_t i = upper_bound(covPos.begin(), covPos.end(), x) - covPos.begin() - 1;
        return covCum[i] + (int64_t)covDepth[i] * (x - covPos[i]);
    };
    vector<int64_t> groupAt{cs}, groupCols{0}, groupBases{0};   /* prefix sums over boundaries */
    for (auto &gb : groups) {
        groupAt.push_back(gb.first);
        groupCols.push_back(groupCols.back() + gb.second.width);
        groupBases.push_back(groupBases.back() + gb.second.bases);
    }
    auto groupsUpTo = [&](int64_t x, const vector<int64_t> &sum) {   /* groups at boundaries <= x */
        return sum[upper_bound(groupAt.begin(), groupAt.end(), x) - groupAt.begin() - 1];
    };
    auto countLE = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(upper_bound(v.begin(), v.end(), x) - v.begin()); };
    auto countLT = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(lower_bound(v.begin(), v.end(), x) - v.begin()); };
    auto rowsIn = [&](int64_t s, int64_t e) { return countLT(starts, e) - countLE(ends, s); };   /* chains overlapping [s, e) */

    /* ---- blocks.  A block [s, e) holds the reference columns s..e-1 and the insertion
     * groups at boundaries (s, e]: a group where the block is cut ends the left block. */
    size_t nextChain = 0;
    vector<size_t> active;
    string text, dna;
    auto emit = [&](int64_t s, int64_t e) {
        while (nextChain < chains.size() && chains[nextChain].start() < e) active.push_back(nextChain++);
        vector<size_t> rows;
        size_t w = 0;
        for (size_t k : active) {
            if (chains[k].end() <= s) continue;
            active[w++] = k;
            rows.push_back(k);
        }
        active.resize(w);
        /* column layout */
        vector<int64_t> refCol((size_t)(e - s));
        map<int64_t, int64_t> groupCol;   /* boundary -> first column */
        int64_t ncol = 0;
        auto gIt = groups.upper_bound(s);
        for (int64_t x = s; x < e; x++) {
            if (x > s && gIt != groups.end() && gIt->first == x) {
                groupCol[x] = ncol;
                ncol += gIt->second.width;
                _stats.insertionColumns += gIt->second.width;
                _stats.insertedBases += gIt->second.bases;
                ++gIt;
            }
            refCol[(size_t)(x - s)] = ncol++;
        }
        if (gIt != groups.end() && gIt->first == e) {
            groupCol[e] = ncol;
            ncol += gIt->second.width;
            _stats.insertionColumns += gIt->second.width;
            _stats.insertedBases += gIt->second.bases;
        }
        struct Row {
            const Chain *c;
            int64_t o1, o2;   /* oriented bases [o1, o2) */
        };
        vector<Row> out;
        for (size_t k : rows) {
            const Chain &c = chains[k];
            int64_t o1 = -1, o2 = -1;
            auto span = [&](int64_t a, int64_t b) {
                if (o1 < 0 || a < o1) o1 = a;
                if (b > o2) o2 = b;
            };
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                if (a < b) span(c.ostart[i] + (a - p.srcStart), c.ostart[i] + (b - p.srcStart));
                if (i < c.at.size() && c.at[i] > s && c.at[i] <= e) span(c.oend(i), c.ostart[i + 1]);
            }
            if (o1 >= 0) out.push_back(Row{&c, o1, o2});
        }
        sort(out.begin(), out.end(), [](const Row &x, const Row &y) {
            int n = x.c->genome->name().compare(y.c->genome->name());
            if (n) return n < 0;
            return make_tuple(x.c->seq, x.o1, x.c->reversed) < make_tuple(y.c->seq, y.o1, y.c->reversed);
        });
        string block = "a\n";
        auto sLine = [&](const Genome *g, uint32_t seq, int64_t start, int64_t size, bool rev, const string &t) {
            block += "s\t";
            if (!_opts.onlySequenceNames) block += g->name() + ".";
            block += g->sequenceName(seq);
            block += "\t" + to_string(start) + "\t" + to_string(size) + "\t" + (rev ? "-" : "+") + "\t" +
                     to_string(g->sequenceLength(seq)) + "\t" + t + "\n";
        };
        fetch(_ref, refSeq, s, e, false, dna);
        text.assign((size_t)ncol, '-');
        for (int64_t x = s; x < e; x++) text[(size_t)refCol[(size_t)(x - s)]] = dna[(size_t)(x - s)];
        sLine(_ref, refSeq, s, e - s, false, text);
        for (const Row &r : out) {
            const Chain &c = *r.c;
            fetch(c.genome, c.seq, r.o1, r.o2, c.reversed, dna);
            text.assign((size_t)ncol, '-');
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                for (int64_t x = a; x < b; x++)
                    text[(size_t)refCol[(size_t)(x - s)]] = dna[(size_t)(c.ostart[i] + (x - p.srcStart) - r.o1)];
                if (i < c.at.size() && c.at[i] > s && c.at[i] <= e) {
                    int64_t base = groupCol[c.at[i]], o = c.oend(i);
                    const vector<int> &cols = c.cols[i];
                    for (size_t j = 0; j < cols.size(); j++) text[(size_t)(base + cols[j])] = dna[(size_t)(o + (int64_t)j - r.o1)];
                }
            }
            sLine(c.genome, c.seq, r.o1, r.o2 - r.o1, c.reversed, text);
            _stats.alignedBases += r.o2 - r.o1;
        }
        block += "\n";
        _out.write(block);
        _stats.blocks++;
        _stats.rows += out.size() + 1;
        _stats.columns += ncol;
    };
    /* close [s, e), splitting at the caps (counted in reference columns) */
    auto close = [&](int64_t s, int64_t e) {
        while (s < e) {
            int64_t w = min(e - s, _opts.maxBlockCols);
            while (w > 1 && (rowsIn(s, s + w) + 1) * w > _opts.maxBlockCells) w = max<int64_t>(1, w / 2);
            emit(s, s + w);
            s += w;
        }
    };

    vector<int64_t> bps;
    for (int64_t x : starts)
        if (x > cs && x < ce) bps.push_back(x);
    for (int64_t x : ends)
        if (x > cs && x < ce) bps.push_back(x);
    sort(bps.begin(), bps.end());
    bps.erase(unique(bps.begin(), bps.end()), bps.end());
    int64_t s = cs;
    for (size_t i = 0; i < bps.size(); i++) {
        int64_t b = bps[i], nb = i + 1 < bps.size() ? bps[i + 1] : ce;
        int64_t endsAt = countLE(ends, b) - countLT(ends, b), startsAt = countLE(starts, b) - countLT(starts, b);
        int64_t activeRows = (countLT(starts, b) - countLT(ends, b)) + startsAt;
        bool cut = activeRows > 0 && (double)(endsAt + startsAt) >= _opts.breakFrac * (double)activeRows;
        if (!cut) {
            /* the block continued to the next breakpoint: its cells and aligned bases,
             * insertion groups included */
            int64_t width = (nb - s) + groupsUpTo(nb, groupCols) - groupsUpTo(s, groupCols);
            int64_t nongap = aligned(nb) - aligned(s) + groupsUpTo(nb, groupBases) - groupsUpTo(s, groupBases);
            int64_t gaps = rowsIn(s, nb) * width - nongap;
            cut = (double)gaps > _opts.gapBudget * (double)nongap;
        }
        if (cut) {
            close(s, b);
            s = b;
        }
    }
    close(s, ce);
}
}
