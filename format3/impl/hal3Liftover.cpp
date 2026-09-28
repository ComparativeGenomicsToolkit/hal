/*
 * HAL format 3: liftover.  See hal3Liftover.h (and hal3BedLiftover.cpp for writing
 * BED as halLiftover does).
 */
#include "hal3Liftover.h"

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>

using namespace std;

namespace hal3 {

/* ---------------------------------------------------------------- mapping */

namespace {

/* The sub-piece of c covering [lo, lo+len) of the current genome, carried to the next
 * genome at [nstart, nstart+len) of sequence nseq, with relative strand s. */
LiftPiece carry(const LiftPiece &c, int64_t lo, int64_t len, uint32_t nseq, int64_t nstart, bool s) {
    int64_t off = c.reversed ? (c.tgtStart + c.length) - (lo + len) : lo - c.tgtStart;
    return LiftPiece{c.srcStart + off, nseq, nstart, len, c.reversed != s};
}

/* The copies of a parent interval [lo, hi) in child -- hits, the axis-b query result --
 * reduced to one per parent base: the copy with the smallest child genome coordinate,
 * which is HAL's canonical paralog in Cactus alignments.  emit(plo, phi, hit) for each
 * stretch of the parent chosen from one hit.  Elementary intervals between hit
 * boundaries are resolved by the concavity of the minimum of lines of slope +-1: if
 * one hit is the minimum at both ends it is throughout, otherwise split. */
void selectCanonical(const Genome *child, const vector<lmap_hit> &hits,
                     const function<void(int64_t, int64_t, const lmap_hit &)> &emit) {
    size_t n = hits.size();
    if (n == 0) return;
    if (n == 1) {
        emit(hits[0].b, hits[0].b + hits[0].len, hits[0]);
        return;
    }
    auto childCoord = [&](const lmap_hit &h, int64_t p) {
        int64_t off = p - h.b;
        return child->sequenceStart(h.a_member) + (h.strand ? h.a + h.len - 1 - off : h.a + off);
    };
    vector<int64_t> cuts;
    for (const lmap_hit &h : hits) {
        cuts.push_back(h.b);
        cuts.push_back(h.b + h.len);
    }
    sort(cuts.begin(), cuts.end());
    cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());
    auto best = [&](int64_t p) {
        int bi = -1;
        int64_t bv = 0;
        for (size_t k = 0; k < n; k++) {
            const lmap_hit &h = hits[k];
            if (p < h.b || p >= h.b + h.len) continue;
            int64_t v = childCoord(h, p);
            if (bi < 0 || v < bv) {
                bi = (int)k;
                bv = v;
            }
        }
        return bi;
    };
    function<void(int64_t, int64_t)> resolve = [&](int64_t lo, int64_t hi) {
        int k0 = best(lo), k1 = best(hi - 1);
        if (k0 < 0) return;
        if (k0 != k1) {
            int64_t mid = lo + (hi - lo) / 2;
            resolve(lo, mid);
            resolve(mid, hi);
            return;
        }
        emit(lo, hi, hits[k0]);
    };
    for (size_t k = 0; k + 1 < cuts.size(); k++) resolve(cuts[k], cuts[k + 1]);
}

/* Do two clipped hits come from the same run (same child member, strand and diagonal)? */
bool sameRun(const lmap_hit &x, const lmap_hit &y) {
    if (x.a_member != y.a_member || x.b_member != y.b_member || x.strand != y.strand) return false;
    return x.strand ? x.a + x.b + x.len == y.a + y.b + y.len : x.a - x.b == y.a - y.b;
}
}

void liftUpBranch(const Genome *g, const vector<LiftPiece> &in, vector<LiftPiece> &out, vector<LiftPiece> *unmapped,
                  CursorCache &cursors, bool canonicalOnly) {
    out.clear();
    if (unmapped) unmapped->clear();
    vector<lmap_hit> hits, copies;
    for (const LiftPiece &c : in) {
        cursors.query(g, 0, c.tgtSeq, c.tgtStart, c.tgtStart + c.length, hits);   /* sorted on a */
        int64_t pos = c.tgtStart;
        auto keepUnmapped = [&](int64_t lo, int64_t hi) {
            if (unmapped && hi > lo) unmapped->push_back(carry(c, lo, hi - lo, c.tgtSeq, lo, false));
        };
        for (const lmap_hit &h : hits) {
            keepUnmapped(pos, h.a);
            pos = h.a + h.len;
            if (!canonicalOnly) {
                out.push_back(carry(c, h.a, h.len, h.b_member, h.b, h.strand != 0));
                continue;
            }
            /* only where this base is the canonical copy of its parent base; elsewhere
             * the base stays here, as if it had no parent */
            cursors.query(g, 1, h.b_member, h.b, h.b + h.len, copies);
            int64_t upTo = h.a;     /* child coordinates of h handled so far (forward) */
            vector<pair<int64_t, int64_t>> mine;   /* child intervals where h is canonical */
            selectCanonical(g, copies, [&](int64_t plo, int64_t phi, const lmap_hit &k) {
                if (!sameRun(k, h)) return;
                int64_t alo = h.strand ? h.a + (h.b + h.len - phi) : h.a + (plo - h.b);
                mine.push_back({alo, alo + (phi - plo)});
            });
            sort(mine.begin(), mine.end());
            for (auto &m : mine) {
                keepUnmapped(upTo, m.first);
                int64_t blo = h.strand ? h.b + (h.a + h.len - m.second) : h.b + (m.first - h.a);
                out.push_back(carry(c, m.first, m.second - m.first, h.b_member, blo, h.strand != 0));
                upTo = m.second;
            }
            keepUnmapped(upTo, h.a + h.len);
        }
        keepUnmapped(pos, c.tgtStart + c.length);
    }
}

void liftDownBranch(const Genome *child, bool dupes, const vector<LiftPiece> &in, vector<LiftPiece> &out,
                    CursorCache &cursors) {
    out.clear();
    vector<lmap_hit> hits;
    for (const LiftPiece &c : in) {
        cursors.query(child, 1, c.tgtSeq, c.tgtStart, c.tgtStart + c.length, hits);
        if (dupes) {
            for (const lmap_hit &h : hits) out.push_back(carry(c, h.b, h.len, h.a_member, h.a, h.strand != 0));
            continue;
        }
        selectCanonical(child, hits, [&](int64_t lo, int64_t hi, const lmap_hit &h) {
            int64_t off = lo - h.b, len = hi - lo;
            int64_t a = h.strand ? h.a + h.len - off - len : h.a + off;
            out.push_back(carry(c, lo, len, h.a_member, a, h.strand != 0));
        });
    }
}

CursorCache::CursorCache(size_t bytesPerCursor, size_t maxCursors) : _bytes(bytesPerCursor), _max(maxCursors) {}

CursorCache::~CursorCache() {
    for (auto &e : _open) lmap_cursor_close(e.second.cursor);
}

void CursorCache::query(const Genome *g, int axis, uint32_t m, int64_t lo, int64_t hi, vector<lmap_hit> &out) {
    out.clear();
    lmap_file *f = g->hgn();
    Key key = make_tuple(f, axis, m);
    auto it = _open.find(key);
    if (it == _open.end()) {
        if (_open.size() >= _max && !_lru.empty()) {
            auto oldest = _open.find(_lru.back());
            lmap_cursor_close(oldest->second.cursor);
            _open.erase(oldest);
            _lru.pop_back();
        }
        lmap_cursor *c = lmap_cursor_open(f, axis, m, nullptr, 0, _bytes);
        if (!c) throw runtime_error("cannot open a cursor on " + g->name());
        _lru.push_front(key);
        it = _open.emplace(key, Entry{c, _lru.begin()}).first;
    } else if (it->second.use != _lru.begin()) {
        _lru.splice(_lru.begin(), _lru, it->second.use);
    }
    lmap_cursor *c = it->second.cursor;
    lmap_hit h, clipped;
    int rc = lmap_cursor_seek(c, lo, hi);
    while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
        rc = 0;
        if (lmap_hit_clip(&h, axis, lo, hi, &clipped)) out.push_back(clipped);
    }
    if (rc != 0) throw runtime_error("query failed on " + g->name());
}

void liftInterval(const Genome *src, uint32_t srcSeq, int64_t start, int64_t end, const Genome *tgt, bool dupes,
                  vector<LiftPiece> &out, CursorCache *cursors) {
    out.clear();
    unique_ptr<CursorCache> own;
    if (!cursors) {
        own.reset(new CursorCache());
        cursors = own.get();
    }
    if (end <= start) return;
    const Genome *mrca = lowestCommonAncestor(src, tgt);
    vector<LiftPiece> cur{LiftPiece{start, srcSeq, start, end - start, false}}, next;
    for (const Genome *g = src; g != mrca && !cur.empty(); g = g->parent()) {
        liftUpBranch(g, cur, next, nullptr, *cursors);
        cur.swap(next);
    }
    vector<const Genome *> down;
    for (const Genome *g = tgt; g != mrca; g = g->parent()) down.push_back(g);
    for (auto it = down.rbegin(); it != down.rend() && !cur.empty(); ++it) {
        liftDownBranch(*it, dupes, cur, next, *cursors);
        cur.swap(next);
    }
    out.swap(cur);
}
}
