/*
 * HAL format 3: the universal column axis.  See hal3Universal.h.
 */
#include "hal3Universal.h"
#include "hal3Walk.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#ifdef __GLIBC__
#include <malloc.h>
#endif

using namespace std;

namespace hal3 {

bool placementLess(const Placement &x, const Placement &y) {
    if (x.slot != y.slot) return x.slot < y.slot;
    if (x.right() != y.right()) return !x.right();
    int gx = x.genome(), gy = y.genome();
    if (gx != gy) return x.right() ? gx < gy : gx > gy;
    return x.right() ? x.jcode < y.jcode : x.jcode > y.jcode;
}

void Moments::add(int64_t a, int64_t b) {
    if (a >= b) return;
    /* sums over 0 .. m-1, exact in 128 bits (m < 2^40 here), then taken mod 2^64 */
    auto sum1 = [](unsigned __int128 m) { return m * (m - 1) / 2; };
    auto sum2 = [](unsigned __int128 m) { return m == 0 ? m : (m - 1) * m * (2 * m - 1) / 6; };
    n += (uint64_t)(b - a);
    s1 += (uint64_t)(sum1((uint64_t)b) - sum1((uint64_t)a));
    s2 += (uint64_t)(sum2((uint64_t)b) - sum2((uint64_t)a));
}

void forEachTopInterval(const Genome *g, uint32_t seq, bool canonicalOnly, CursorCache &cursors,
                        const function<void(int64_t, int64_t)> &emit) {
    const int64_t len = g->sequenceLength(seq);
    int64_t s = 0, e = 0;   /* the interval being grown, [s, e); none while s == e */
    auto add = [&](int64_t a, int64_t b) {
        if (a >= b) return;
        if (a != e) {
            if (e > s) emit(s, e);
            s = a;
        }
        e = b;
    };
    if (!canonicalOnly) {
        /* the gaps between the branch's runs (which never overlap on this side): one
         * streaming pass over them, nothing kept */
        unique_ptr<lmap_cursor, void (*)(lmap_cursor *)> c(lmap_cursor_open(g->hgn(), 0, seq, nullptr, 0, 0),
                                                            lmap_cursor_close);
        if (!c) throw runtime_error("cannot open a cursor on " + g->name());
        lmap_hit h;
        int64_t pos = 0;
        int rc = lmap_cursor_seek(c.get(), 0, len);
        while (rc == 0 && (rc = lmap_cursor_next(c.get(), &h)) == 1) {
            rc = 0;
            add(pos, h.a);
            pos = max(pos, h.a + h.len);
        }
        if (rc != 0) throw runtime_error("cannot read the branch above " + g->name());
        add(pos, len);
    } else {
        /* what the canonical-only lift up leaves behind, a window at a time */
        const int64_t WINDOW = (int64_t)4 << 20;
        vector<LiftPiece> in(1), up, left;
        for (int64_t w = 0; w < len; w += WINDOW) {
            in[0] = LiftPiece{w, seq, w, min(len, w + WINDOW) - w, false};
            liftUpBranch(g, in, up, &left, cursors, true);
            sort(left.begin(), left.end(), [](const LiftPiece &x, const LiftPiece &y) { return x.tgtStart < y.tgtStart; });
            for (const LiftPiece &p : left) add(p.tgtStart, p.tgtStart + p.length);
        }
    }
    if (e > s) emit(s, e);
}

/* ---------------------------------------------------------------- the index */

namespace {

/* Placements appended in blocks that double up to 64k entries: no vector's slack of up
 * to 2x while the index is built, and freed block by block as it is gathered. */
class Blocks {
  public:
    void push(const Placement &p) {
        if (_blocks.empty() || _blocks.back().size() == _blocks.back().capacity()) {
            size_t cap = _blocks.empty() ? 64 : min<size_t>(2 * _blocks.back().capacity(), (size_t)1 << 16);
            _blocks.emplace_back();
            _blocks.back().reserve(cap);
        }
        _blocks.back().push_back(p);
        _size++;
    }
    size_t size() const { return _size; }
    /* Append everything to v, emptying this.  Freed blocks stay with the heap, which a
     * large v does not reuse: their pages go back every 64 MB (moved counts), or the
     * index would take twice its size while it is gathered. */
    void moveTo(vector<Placement> &v, size_t &moved) {
        for (vector<Placement> &b : _blocks) {
            v.insert(v.end(), b.begin(), b.end());
            moved += b.size() * sizeof(Placement);
            vector<Placement>().swap(b);
#ifdef __GLIBC__
            if (moved >= ((size_t)64 << 20)) {
                malloc_trim(0);
                moved = 0;
            }
#endif
        }
        _blocks.clear();
        _size = 0;
    }

  private:
    vector<vector<Placement>> _blocks;
    size_t _size = 0;
};

/* sort v on `threads` threads: parts sorted at once, then merged */
void parallelSort(vector<Placement> &v, int threads) {
    size_t parts = v.size() < ((size_t)1 << 20) ? 1 : (size_t)max(1, threads);
    if (parts == 1) {
        sort(v.begin(), v.end(), placementLess);
        return;
    }
    size_t step = (v.size() + parts - 1) / parts;
    orderedParallel<char>(parts, threads,
                          [&](size_t i) {
                              size_t lo = min(v.size(), i * step), hi = min(v.size(), lo + step);
                              sort(v.begin() + (ptrdiff_t)lo, v.begin() + (ptrdiff_t)hi, placementLess);
                              return (char)0;
                          },
                          [](size_t, char &) {});
    for (size_t w = step; w < v.size(); w *= 2)
        for (size_t lo = 0; lo + w < v.size(); lo += 2 * w)
            inplace_merge(v.begin() + (ptrdiff_t)lo, v.begin() + (ptrdiff_t)(lo + w),
                          v.begin() + (ptrdiff_t)min(v.size(), lo + 2 * w), placementLess);
}
}

UniversalIndex::UniversalIndex(const Alignment &aln, const Genome *backbone, bool leafInsertions, bool dupes,
                               const set<const Genome *> &excluded, bool allTops, int threads) {
    /* the scope, in the Newick's preorder (any genome it does not name after, as found) */
    set<const Genome *> scope;
    vector<const Genome *> found;
    for (vector<const Genome *> st{backbone}; !st.empty();) {
        const Genome *g = st.back();
        st.pop_back();
        scope.insert(g);
        found.push_back(g);
        for (auto c = g->children().rbegin(); c != g->children().rend(); ++c) st.push_back(*c);
    }
    auto enroll = [&](const Genome *g) {
        if (!scope.count(g) || _rank.count(g)) return;
        _rank[g] = (int)_genomes.size();
        _genomes.push_back(g);
    };
    for (const NewickNode &n : parseNewick(aln.manifest().newick)) {
        const Genome *g = aln.genome(n.name);
        if (g) enroll(g);
    }
    for (const Genome *g : found) enroll(g);
    if (_genomes[0] != backbone) throw runtime_error("internal error: the backbone is not first in preorder");
    const size_t ng = _genomes.size();
    if (ng >= ((size_t)1 << 22)) throw runtime_error("too many genomes for --universal");
    _covered.assign(ng, 0);
    _codeStart.resize(ng);
    _anchored.resize(ng);
    _topIntervals.assign(ng, 0);
    _topBases.resize(ng);
    _allBases.resize(ng);
    for (size_t r = 0; r < ng; r++) {
        const Genome *g = _genomes[r];
        _covered[r] = r > 0 && !excluded.count(g) && (leafInsertions || !g->children().empty());
        for (uint32_t q = 0; q < g->numSequences(); q++) _codeStart[r].push_back(g->sequenceStart(q) + q);
        int64_t total = g->numSequences() ? g->sequenceStart(g->numSequences() - 1) + g->sequenceLength(g->numSequences() - 1) : 0;
        _allBases[r].add(0, total);
    }
    _topBases[0] = _allBases[0];
    _topIntervals[0] = backbone->numSequences();

    /* tasks: runs of whole sequences of about 20 Mb */
    struct Task {
        int g;
        uint32_t q0, q1;
    };
    vector<Task> tasks;
    for (size_t r = 1; r < ng; r++) {
        if (!_covered[r] && !(allTops && !excluded.count(_genomes[r]))) continue;
        const Genome *g = _genomes[r];
        int64_t acc = 0;
        uint32_t q0 = 0;
        for (uint32_t q = 0; q < g->numSequences(); q++) {
            acc += g->sequenceLength(q);
            if (acc >= 20000000 || q + 1 == g->numSequences()) {
                tasks.push_back(Task{(int)r, q0, q + 1});
                q0 = q + 1;
                acc = 0;
            }
        }
    }
    const int nt = max(1, threads);
    struct Worker {
        unique_ptr<CursorCache> cursors;
        vector<Blocks> anchored;
        vector<Placement> flankless;
        vector<uint64_t> intervals;
        vector<Moments> bases;
    };
    vector<Worker> workers((size_t)nt);
    orderedParallel<char>(
        tasks.size(), nt,
        [&](size_t k, int w) {
            Worker &W = workers[(size_t)w];
            if (!W.cursors) {
                W.cursors.reset(new CursorCache((size_t)1 << 20, 256));
                W.anchored.resize(ng);
                W.intervals.assign(ng, 0);
                W.bases.resize(ng);
            }
            const Task &t = tasks[k];
            const Genome *G = _genomes[(size_t)t.g];
            const bool cov = _covered[(size_t)t.g] != 0;
            struct J {
                int64_t jcode, len;
                bool rightFlank;   /* placed by R (before it), not by L (after it) */
            };
            vector<J> js;
            vector<LiftPiece> cur, next, left;
            for (uint32_t q = t.q0; q < t.q1; q++) {
                const int64_t len = G->sequenceLength(q), gs = G->sequenceStart(q);
                forEachTopInterval(G, q, !dupes, *W.cursors, [&](int64_t s, int64_t e) {
                    W.intervals[(size_t)t.g]++;
                    W.bases[(size_t)t.g].add(gs + s, gs + e);
                    if (!cov) return;
                    if (s > 0 || e < len) {
                        cur.push_back(LiftPiece{(int64_t)js.size(), q, s > 0 ? s - 1 : e, 1, false});
                        js.push_back(J{code(t.g, q, s), e - s, s == 0});
                    } else {
                        W.flankless.push_back(Placement{-1, code(t.g, q, 0), Placement::pack(len, t.g, false, false)});
                    }
                });
            }
            /* each flank up to its top base K, in A: the slot next to K, as the table in
             * hal3Universal.h's order describes -- after L in G's direction, or before R */
            for (const Genome *a = G; !cur.empty(); a = a->parent()) {
                const int ar = rank(a);
                auto place = [&](const LiftPiece &p) {
                    const J &j = js[(size_t)p.srcStart];
                    const bool rev = p.reversed;
                    const int64_t k = j.rightFlank ? (rev ? p.tgtStart + 1 : p.tgtStart) : (rev ? p.tgtStart : p.tgtStart + 1);
                    W.anchored[(size_t)ar].push(
                        Placement{code(ar, p.tgtSeq, k), j.jcode, Placement::pack(j.len, t.g, rev, rev != j.rightFlank)});
                };
                if (a == backbone) {
                    for (const LiftPiece &p : cur) place(p);
                    break;
                }
                sort(cur.begin(), cur.end(), [](const LiftPiece &x, const LiftPiece &y) {
                    return x.tgtSeq != y.tgtSeq ? x.tgtSeq < y.tgtSeq : x.tgtStart < y.tgtStart;
                });
                liftUpBranch(a, cur, next, &left, *W.cursors, !dupes);
                if (a == G && !left.empty())
                    throw runtime_error("internal error: a flank of a top interval of " + G->name() + " has no parent");
                for (const LiftPiece &p : left) place(p);
                cur.swap(next);
            }
            return (char)0;
        },
        [](size_t, char &) {});

    /* gather, freeing the blocks as it goes, and sort */
    size_t moved = 0;
    for (size_t r = 0; r < ng; r++) {
        size_t n = 0;
        for (Worker &W : workers)
            if (W.cursors) n += W.anchored[r].size();
        vector<Placement> &v = _anchored[r];
        v.reserve(n);
        for (Worker &W : workers)
            if (W.cursors) W.anchored[r].moveTo(v, moved);
        parallelSort(v, nt);
    }
    for (Worker &W : workers) {
        if (!W.cursors) continue;
        _flankless.insert(_flankless.end(), W.flankless.begin(), W.flankless.end());
        for (size_t r = 1; r < ng; r++) {
            _topIntervals[r] += W.intervals[r];
            _topBases[r].add(W.bases[r]);
        }
    }
    workers.clear();
    sort(_flankless.begin(), _flankless.end(), [](const Placement &x, const Placement &y) {
        return x.genome() != y.genome() ? x.genome() < y.genome() : x.jcode < y.jcode;
    });

    /* every covered interval placed once: counts and bases agree */
    uint64_t want = 0, bases = 0, placedBases = 0;
    for (size_t r = 1; r < ng; r++)
        if (_covered[r]) {
            want += _topIntervals[r];
            bases += _topBases[r].n;
        }
    _placements = _flankless.size();
    for (const Placement &p : _flankless) placedBases += (uint64_t)p.len();
    for (size_t r = 0; r < ng; r++) {
        _placements += _anchored[r].size();
        for (const Placement &p : _anchored[r]) placedBases += (uint64_t)p.len();
        if (!_anchored[r].empty() && _genomes[r]->children().empty())
            throw runtime_error("internal error: an interval anchored in leaf " + _genomes[r]->name());
    }
    if (_placements != want || placedBases != bases)
        throw runtime_error("internal error: " + to_string(_placements) + " intervals placed (" + to_string(placedBases) +
                            " bases) of " + to_string(want) + " (" + to_string(bases) + ")");
    _columns = (int64_t)(_allBases[0].n + bases);
}

int UniversalIndex::rank(const Genome *g) const {
    auto r = _rank.find(g);
    return r == _rank.end() ? -1 : r->second;
}

void UniversalIndex::decode(int rank, int64_t code, uint32_t &q, int64_t &k) const {
    const vector<int64_t> &cs = _codeStart[(size_t)rank];
    q = (uint32_t)(upper_bound(cs.begin(), cs.end(), code) - cs.begin() - 1);
    k = code - cs[q];
}

/* ---------------------------------------------------------------- the planner */

UniversalPlanner::UniversalPlanner(const UniversalIndex &index, const vector<TopRange> &backbone, bool whole, int64_t chunk)
    : _idx(index), _ranges(backbone), _whole(whole), _chunk(chunk) {
    if (whole)
        for (const Placement &p : index.flankless()) {
            uint32_t q;
            int64_t k;
            index.decode(p.genome(), p.jcode, q, k);
            _ranges.push_back(TopRange{p.genome(), q, 0, p.len()});
        }
}

/* The run [lo, hi) of genome g, with the placements anchored on its bases: slots lo .. hi,
 * less a slot-lo placement anchored on lo-1 and a slot-hi one anchored on hi. */
UniversalPlanner::Frame UniversalPlanner::frame(int g, int64_t lo, int64_t hi, bool rev, int depth) const {
    const vector<Placement> &a = _idx.anchoredIn(g);
    size_t b = (size_t)(lower_bound(a.begin(), a.end(), lo, [](const Placement &p, int64_t x) { return p.slot < x; }) - a.begin());
    size_t e = (size_t)(upper_bound(a.begin(), a.end(), hi, [](int64_t x, const Placement &p) { return x < p.slot; }) - a.begin());
    while (b < e && a[b].slot == lo && !a[b].right()) b++;
    while (e > b && a[e - 1].slot == hi && a[e - 1].right()) e--;
    return Frame{g, depth, rev, lo, hi, rev ? hi : lo, b, e, rev ? e : b};
}

/* Up to `want` more columns of the walk on stack st (to out, if given); returns how many.
 * A placement is entered as a frame of its own (visited counts them), laid out reversed
 * when its strand differs from the orientation of the run it hangs off. */
int64_t UniversalPlanner::step(vector<Frame> &st, int64_t want, vector<Segment> *out, uint64_t &visited, int &maxDepth) const {
    int64_t cols = 0;
    while (cols < want && !st.empty()) {
        Frame &f = st.back();
        const vector<Placement> &a = _idx.anchoredIn(f.g);
        const Placement *p = nullptr;
        if (!f.rev) {
            if (f.it < f.e && a[f.it].slot == f.cur) p = &a[f.it++];
        } else if (f.it > f.b && a[f.it - 1].slot == f.cur) {
            p = &a[--f.it];
        }
        if (p) {
            const int d = f.depth + 1;
            const bool r = p->rev() != f.rev;
            st.push_back(frame(p->genome(), p->jcode, p->jcode + p->len(), r, d));   /* f is invalid from here */
            visited++;
            maxDepth = max(maxDepth, d);
            continue;
        }
        if (f.cur == (f.rev ? f.lo : f.hi)) {
            st.pop_back();
            continue;
        }
        const int64_t next = !f.rev ? (f.it < f.e ? a[f.it].slot : f.hi) : (f.it > f.b ? a[f.it - 1].slot : f.lo);
        const int64_t k = min(f.rev ? f.cur - next : next - f.cur, want - cols);
        if (out) {
            uint32_t q;
            int64_t pos;
            _idx.decode(f.g, f.rev ? f.cur - k : f.cur, q, pos);
            out->push_back(Segment{cols, f.g, q, pos, k, f.rev, f.depth});
        }
        f.cur += f.rev ? -k : k;
        cols += k;
    }
    return cols;
}

void UniversalPlanner::run() {
    _chunks.clear();
    _snaps.clear();
    _u = 0;
    _visited = 0;
    _maxDepth = 0;
    const size_t nb = _ranges.size() - (_whole ? _idx.flankless().size() : 0);
    for (size_t r = 0; r < _ranges.size(); r++) {
        const TopRange &t = _ranges[r];
        if (r >= nb) _visited++;   /* a flankless interval */
        vector<Frame> st{frame(t.g, _idx.code(t.g, t.seq, t.start), _idx.code(t.g, t.seq, t.end), false, 0)};
        for (;;) {
            Chunk c{r, _u, 0, t.start + (st.empty() ? 0 : st[0].cur - st[0].lo)};
            vector<Frame> snap = st;
            c.cols = step(st, _chunk, nullptr, _visited, _maxDepth);
            if (c.cols == 0) break;
            _chunks.push_back(c);
            _snaps.push_back(move(snap));
            _u += c.cols;
        }
    }
    if (_whole && (_visited != _idx.placements() || _u != _idx.columns()))
        throw runtime_error("internal error: the universal axis visited " + to_string(_visited) + " of " +
                            to_string(_idx.placements()) + " intervals, " + to_string(_u) + " of " +
                            to_string(_idx.columns()) + " columns");
}

void UniversalPlanner::segments(size_t k, vector<Segment> &out) const {
    out.clear();
    vector<Frame> st = _snaps[k];
    uint64_t visited = 0;
    int depth = 0;
    if (step(st, _chunks[k].cols, &out, visited, depth) != _chunks[k].cols)
        throw runtime_error("internal error: chunk " + to_string(k) + " of the universal axis replays short");
}
}
