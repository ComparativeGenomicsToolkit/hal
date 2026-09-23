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

MafExporter::MafExporter(const Alignment &alignment, const Genome *ref, const MafOptions &opts, MafSink &out)
    : _aln(alignment), _ref(ref), _opts(opts), _out(out), _cursors((size_t)1 << 20, 4096) {
    /* printed genomes: the targets (or the subtree of rootGenome) and the reference;
     * scope: the spanning tree of those, whose top is where the walk turns */
    vector<const Genome *> targets;
    if (!_opts.rootGenome.empty()) {
        if (!_opts.targets.empty()) throw runtime_error("--rootGenome and --targetGenomes cannot be combined");
        const Genome *r = _aln.genome(_opts.rootGenome);
        if (!r) throw runtime_error("no genome " + _opts.rootGenome);
        vector<const Genome *> stack{r};
        while (!stack.empty()) {
            const Genome *g = stack.back();
            stack.pop_back();
            targets.push_back(g);
            for (const Genome *c : g->children()) stack.push_back(c);
        }
    }
    for (const string &t : _opts.targets) {
        const Genome *g = _aln.genome(t);
        if (!g) throw runtime_error("no genome " + t);
        targets.push_back(g);
    }
    if (targets.empty()) {
        vector<const Genome *> stack{_aln.root()};
        while (!stack.empty()) {
            const Genome *g = stack.back();
            stack.pop_back();
            targets.push_back(g);
            for (const Genome *c : g->children()) stack.push_back(c);
        }
    }
    targets.push_back(_ref);
    _top = targets[0];
    for (const Genome *t : targets) _top = lowestCommonAncestor(_top, t);
    for (const Genome *t : targets) {
        for (const Genome *g = t; g != _top; g = g->parent()) _scope.insert(g);
        _printed.insert(t);
    }
    _scope.insert(_top);
    if (!_opts.ancestors) {
        if (!_ref->children().empty()) throw runtime_error("--noAncestors needs a leaf reference genome");
        for (auto it = _printed.begin(); it != _printed.end();)
            it = (*it)->children().empty() ? next(it) : _printed.erase(it);
    }
}

MafExporter::~MafExporter() {
    for (auto &d : _dna)
        if (d.second) tb_close(d.second);
}

void MafExporter::writeHeader() { _out.write("##maf version=1 scoring=N/A\n# hal " + _aln.manifest().newick + "\n\n"); }

void MafExporter::exportInterval(uint32_t refSeq, int64_t start, int64_t end) {
    for (int64_t s = start; s < end; s += _opts.chunk) chunk(refSeq, s, min(end, s + _opts.chunk));
}

/* ---------------------------------------------------------------- the walk */

namespace {

/* Merge pieces contiguous in both the reference and the target genome. */
void coalesce(vector<LiftPiece> &v) {
    if (v.size() < 2) return;
    auto diag = [](const LiftPiece &p) { return p.reversed ? p.srcStart + p.tgtStart + p.length : p.tgtStart - p.srcStart; };
    sort(v.begin(), v.end(), [&](const LiftPiece &a, const LiftPiece &b) {
        return make_tuple(a.tgtSeq, a.reversed, diag(a), a.srcStart) < make_tuple(b.tgtSeq, b.reversed, diag(b), b.srcStart);
    });
    vector<LiftPiece> o;
    for (const LiftPiece &p : v) {
        if (!o.empty()) {
            LiftPiece &q = o.back();
            if (q.tgtSeq == p.tgtSeq && q.reversed == p.reversed && diag(q) == diag(p) && q.srcStart + q.length == p.srcStart) {
                if (q.reversed) q.tgtStart = p.tgtStart;
                q.length += p.length;
                continue;
            }
        }
        o.push_back(p);
    }
    v.swap(o);
}
}

/* At v: its own pieces (arrived from the parent, plus any tops at v), recorded if v is
 * printed, then carried into each child in scope but skipChild. */
void MafExporter::descend(const Genome *v, vector<LiftPiece> &in, const Genome *skipChild,
                          map<const Genome *, vector<LiftPiece>> &tops, map<const Genome *, vector<LiftPiece>> &rows) {
    auto t = tops.find(v);
    if (t != tops.end()) {
        in.insert(in.end(), t->second.begin(), t->second.end());
        tops.erase(t);
    }
    coalesce(in);
    if (_printed.count(v)) {
        vector<LiftPiece> &r = rows[v];
        r.insert(r.end(), in.begin(), in.end());
    }
    for (const Genome *c : v->children()) {
        if (c == skipChild || !_scope.count(c)) continue;
        vector<LiftPiece> down;
        liftDownBranch(c, _opts.dupes, in, down, _cursors);
        if (!down.empty() || tops.count(c)) descend(c, down, nullptr, tops, rows);
    }
}

void MafExporter::walk(uint32_t refSeq, int64_t start, int64_t end, map<const Genome *, vector<LiftPiece>> &rows) {
    rows.clear();
    vector<LiftPiece> cur{LiftPiece{start, refSeq, start, end - start, false}}, next, unmapped;
    if (!_opts.onlyOrthologs) {
        /* up to the scope top; every top descends into its whole subtree */
        map<const Genome *, vector<LiftPiece>> tops;
        const Genome *a = _ref;
        for (;;) {
            if (a == _top) {
                vector<LiftPiece> &t = tops[a];
                t.insert(t.end(), cur.begin(), cur.end());
                break;
            }
            liftUpBranch(a, cur, next, &unmapped, _cursors, !_opts.dupes);
            if (!unmapped.empty()) {
                vector<LiftPiece> &t = tops[a];
                t.insert(t.end(), unmapped.begin(), unmapped.end());
            }
            cur.swap(next);
            if (cur.empty()) break;
            a = a->parent();
        }
        vector<LiftPiece> none;
        descend(_top, none, nullptr, tops, rows);
    } else {
        /* liftover to the common ancestor with each genome: at each ancestor on the
         * reference's lineage, descend only into the children off the lineage; and
         * into the reference's own subtree */
        map<const Genome *, vector<LiftPiece>> noTops;
        const Genome *below = nullptr;
        for (const Genome *a = _ref;; a = a->parent()) {
            if (a != _ref && _printed.count(a)) {
                vector<LiftPiece> &r = rows[a];
                r.insert(r.end(), cur.begin(), cur.end());
            }
            for (const Genome *c : a->children()) {
                if (c == below || !_scope.count(c)) continue;
                vector<LiftPiece> down;
                liftDownBranch(c, _opts.dupes, cur, down, _cursors);
                if (!down.empty()) descend(c, down, nullptr, noTops, rows);
            }
            if (a == _top) break;
            liftUpBranch(a, cur, next, nullptr, _cursors, !_opts.dupes);
            cur.swap(next);
            coalesce(cur);
            if (cur.empty()) break;
            below = a;
        }
    }
    /* the reference's own identity pieces: its anchor row stands for them */
    auto r = rows.find(_ref);
    if (r != rows.end()) {
        vector<LiftPiece> &v = r->second;
        v.erase(remove_if(v.begin(), v.end(),
                          [&](const LiftPiece &p) { return !p.reversed && p.tgtSeq == refSeq && p.tgtStart == p.srcStart; }),
                v.end());
    }
}

/* ---------------------------------------------------------------- rows and blocks */

/* A row: pieces of one genome, contiguous in it, in reference order.  Oriented
 * coordinates (MAF's: position on the row's strand) increase along the reference. */
struct MafExporter::Chain {
    const Genome *genome;
    uint32_t seq;
    bool reversed;
    vector<LiftPiece> pieces;   /* in reference order */
    vector<int64_t> ostart;     /* oriented start of each piece */
    int64_t start() const { return pieces.front().srcStart; }
    int64_t end() const { return pieces.back().srcStart + pieces.back().length; }
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
    walk(refSeq, cs, ce, lifted);

    /* ---- rows: chain each genome's pieces where they follow each other directly */
    vector<Chain> chains;
    for (auto &gp : lifted) {
        const Genome *g = gp.first;
        vector<LiftPiece> &ps = gp.second;
        size_t n = ps.size();
        vector<int64_t> os(n);
        for (size_t i = 0; i < n; i++) {
            const LiftPiece &p = ps[i];
            os[i] = p.reversed ? g->sequenceLength(p.tgtSeq) - (p.tgtStart + p.length) : p.tgtStart;
        }
        /* candidate links p -> q: q starts where p ends (oriented), later on the reference */
        map<tuple<uint32_t, bool, int64_t>, vector<size_t>> startsAt;
        for (size_t i = 0; i < n; i++) startsAt[make_tuple(ps[i].tgtSeq, ps[i].reversed, os[i])].push_back(i);
        vector<tuple<int64_t, int64_t, size_t, size_t>> edges;   /* dRef, p.src, p, q */
        for (size_t i = 0; i < n; i++) {
            auto it = startsAt.find(make_tuple(ps[i].tgtSeq, ps[i].reversed, os[i] + ps[i].length));
            if (it == startsAt.end()) continue;
            int64_t pend = ps[i].srcStart + ps[i].length;
            for (size_t j : it->second)
                if (ps[j].srcStart >= pend) edges.emplace_back(ps[j].srcStart - pend, ps[i].srcStart, i, j);
        }
        sort(edges.begin(), edges.end());
        vector<int64_t> succ(n, -1), pred(n, -1);
        for (auto &e : edges) {
            size_t p = get<2>(e), q = get<3>(e);
            if (succ[p] < 0 && pred[q] < 0) {
                succ[p] = (int64_t)q;
                pred[q] = (int64_t)p;
            }
        }
        for (size_t i = 0; i < n; i++) {
            if (pred[i] >= 0) continue;
            Chain c{g, ps[i].tgtSeq, ps[i].reversed, {}, {}};
            for (int64_t k = (int64_t)i; k >= 0; k = succ[k]) {
                c.pieces.push_back(ps[k]);
                c.ostart.push_back(os[k]);
            }
            chains.push_back(move(c));
        }
    }
    sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) { return a.start() < b.start(); });

    /* ---- counts for the block rule: chain starts and ends, and covered bases */
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
    vector<int64_t> covPos{cs}, covCum{0};   /* aligned bases before covPos[i] */
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
    auto aligned = [&](int64_t x) {   /* aligned bases in [cs, x) */
        size_t i = upper_bound(covPos.begin(), covPos.end(), x) - covPos.begin() - 1;
        return covCum[i] + (int64_t)covDepth[i] * (x - covPos[i]);
    };
    auto countLE = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(upper_bound(v.begin(), v.end(), x) - v.begin()); };
    auto countLT = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(lower_bound(v.begin(), v.end(), x) - v.begin()); };
    auto rowsIn = [&](int64_t s, int64_t e) { return countLT(starts, e) - countLE(ends, s); };   /* chains overlapping [s, e) */

    /* ---- blocks */
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
        struct Row {
            const Chain *c;
            int64_t o1, o2;   /* oriented bases [o1, o2) */
        };
        vector<Row> out;
        for (size_t k : rows) {
            const Chain &c = chains[k];
            int64_t o1 = -1, o2 = -1;
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                if (a >= b) continue;
                if (o1 < 0) o1 = c.ostart[i] + (a - p.srcStart);
                o2 = c.ostart[i] + (b - p.srcStart);
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
        sLine(_ref, refSeq, s, e - s, false, dna);
        for (const Row &r : out) {
            const Chain &c = *r.c;
            fetch(c.genome, c.seq, r.o1, r.o2, c.reversed, dna);
            text.assign((size_t)(e - s), '-');
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                if (a >= b) continue;
                int64_t off = c.ostart[i] + (a - p.srcStart) - r.o1;
                memcpy(&text[a - s], &dna[off], (size_t)(b - a));
            }
            sLine(c.genome, c.seq, r.o1, r.o2 - r.o1, c.reversed, text);
            _stats.alignedBases += r.o2 - r.o1;
        }
        block += "\n";
        _out.write(block);
        _stats.blocks++;
        _stats.rows += out.size() + 1;
        _stats.columns += e - s;
    };
    /* close [s, e), splitting at the caps */
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
            int64_t width = nb - s, nongap = aligned(nb) - aligned(s);
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
