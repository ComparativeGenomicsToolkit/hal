/*
 * HAL format 3: liftover.  See hal3Liftover.h.
 *
 * The output half follows halLiftover's code closely -- liftover/impl/halLiftover.cpp,
 * halBlockLiftover.cpp and halBlockMapper.cpp (BlockMapper::extractSegment), and
 * api/impl/halMappedSegment.cpp / halSegmentMapper.cpp for the ordering, the merge
 * test and insertAndBreakOverlaps -- because the goal is byte-identical output.
 * Where it differs in form (pieces instead of segment iterators), the comments say
 * which HAL code it stands for.
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

/* A source interval paired with an interval of the genome currently reached. */
struct Cur {
    int64_t srcStart;
    uint32_t seq;
    int64_t start;
    int64_t length;
    bool reversed;
};

/* The sub-piece of c covering [lo, lo+len) of the current genome, carried to the next
 * genome at [nstart, nstart+len) of sequence nseq, with relative strand s. */
Cur carry(const Cur &c, int64_t lo, int64_t len, uint32_t nseq, int64_t nstart, bool s) {
    int64_t off = c.reversed ? (c.start + c.length) - (lo + len) : lo - c.start;
    return Cur{c.srcStart + off, nseq, nstart, len, c.reversed != s};
}

/* Up the branch above g: each base has at most one parent position. */
void stepUp(const Genome *g, vector<Cur> &cur, CursorCache &cc) {
    vector<Cur> next;
    vector<lmap_hit> hits;
    for (const Cur &c : cur) {
        cc.query(g, 0, c.seq, c.start, c.start + c.length, hits);
        for (const lmap_hit &h : hits) next.push_back(carry(c, h.a, h.len, h.b_member, h.b, h.strand != 0));
    }
    cur.swap(next);
}

/* Down the branch above child: every copy, or with dupes off the copy with the
 * smallest child genome coordinate at each parent base. */
void stepDown(const Genome *child, bool dupes, vector<Cur> &cur, CursorCache &cc) {
    vector<Cur> next;
    vector<lmap_hit> hits;
    for (const Cur &c : cur) {
        cc.query(child, 1, c.seq, c.start, c.start + c.length, hits);
        size_t n = hits.size();
        if (dupes || n <= 1) {
            for (const lmap_hit &h : hits) next.push_back(carry(c, h.b, h.len, h.a_member, h.a, h.strand != 0));
            continue;
        }
        const lmap_hit *hg = hits.data();
        /* Elementary intervals of the parent axis between hit boundaries; in each, the
         * candidates' child coordinates are lines of slope +-1, whose minimum is
         * concave -- so if one hit is the minimum at both ends it is throughout, and
         * otherwise the interval is split. */
        auto childCoord = [&](const lmap_hit &h, int64_t p) {
            int64_t off = p - h.b;
            return child->sequenceStart(h.a_member) + (h.strand ? h.a + h.len - 1 - off : h.a + off);
        };
        vector<int64_t> cuts;
        for (size_t k = 0; k < n; k++) {
            cuts.push_back(hg[k].b);
            cuts.push_back(hg[k].b + hg[k].len);
        }
        sort(cuts.begin(), cuts.end());
        cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());
        auto best = [&](int64_t p) {
            int bi = -1;
            int64_t bv = 0;
            for (size_t k = 0; k < n; k++) {
                const lmap_hit &h = hg[k];
                if (p < h.b || p >= h.b + h.len) continue;
                int64_t v = childCoord(h, p);
                if (bi < 0 || v < bv) { bi = (int)k; bv = v; }
            }
            return bi;
        };
        function<void(int64_t, int64_t)> emit = [&](int64_t lo, int64_t hi) {   /* [lo, hi) */
            int k0 = best(lo), k1 = best(hi - 1);
            if (k0 < 0) return;
            if (k0 != k1) {
                int64_t mid = lo + (hi - lo) / 2;
                emit(lo, mid);
                emit(mid, hi);
                return;
            }
            const lmap_hit &h = hg[k0];
            int64_t off = lo - h.b, len = hi - lo;
            int64_t a = h.strand ? h.a + h.len - off - len : h.a + off;
            next.push_back(carry(c, lo, len, h.a_member, a, h.strand != 0));
        };
        for (size_t k = 0; k + 1 < cuts.size(); k++) emit(cuts[k], cuts[k + 1]);
    }
    cur.swap(next);
}
}

CursorCache::CursorCache(size_t bytesPerCursor, size_t maxCursors) : _bytes(bytesPerCursor), _max(maxCursors) {}

CursorCache::~CursorCache() {
    for (auto &e : _open) lmap_cursor_close(e.second.cursor);
}

void CursorCache::query(const Genome *g, int axis, uint32_t m, int64_t lo, int64_t hi, vector<lmap_hit> &out) {
    out.clear();
    lmap_file *f = g->hgn();
    auto key = make_tuple(f, axis, m);
    auto it = _open.find(key);
    if (it == _open.end()) {
        if (_open.size() >= _max) {
            auto oldest = _open.begin();
            for (auto i = _open.begin(); i != _open.end(); ++i)
                if (i->second.lastUse < oldest->second.lastUse) oldest = i;
            lmap_cursor_close(oldest->second.cursor);
            _open.erase(oldest);
        }
        lmap_cursor *c = lmap_cursor_open(f, axis, m, nullptr, 0, _bytes);
        if (!c) throw runtime_error("cannot open a cursor on " + g->name());
        it = _open.emplace(key, Entry{c, 0}).first;
    }
    it->second.lastUse = ++_clock;
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
    vector<Cur> cur{Cur{start, srcSeq, start, end - start, false}};
    for (const Genome *g = src; g != mrca && !cur.empty(); g = g->parent()) stepUp(g, cur, *cursors);
    vector<const Genome *> down;
    for (const Genome *g = tgt; g != mrca; g = g->parent()) down.push_back(g);
    for (auto it = down.rbegin(); it != down.rend() && !cur.empty(); ++it) stepDown(*it, dupes, cur, *cursors);
    for (const Cur &c : cur) out.push_back(LiftPiece{c.srcStart, c.seq, c.start, c.length, c.reversed});
}

/* ---------------------------------------------------------------- halLiftover's output */

namespace {

/* A mapped segment, as halLiftover's MappedSegment: a target interval and a source
 * interval, each with an orientation, in genome coordinates (inclusive ends). */
struct MSeg {
    int64_t tLo, tHi;
    bool tRev;
    int64_t sLo, sHi;
    bool sRev;
    uint32_t tSeq;
    int64_t tStart() const { return tRev ? tHi : tLo; }   /* getStartPosition() */
    int64_t tEnd() const { return tRev ? tLo : tHi; }     /* getEndPosition() */
    int64_t sStart() const { return sRev ? sHi : sLo; }
    int64_t sEnd() const { return sRev ? sLo : sHi; }
    /* source and target run the same way, or opposite ways */
    bool parallel() const { return tRev == sRev; }
};

/* MappedSegment::lessThan (fastComp == slowComp): target (low, high), then source. */
bool msLess(const MSeg &a, const MSeg &b) {
    if (a.tLo != b.tLo) return a.tLo < b.tLo;
    if (a.tHi != b.tHi) return a.tHi < b.tHi;
    if (a.sLo != b.sLo) return a.sLo < b.sLo;
    return a.sHi < b.sHi;
}

/* The part of m on target [lo, hi] (inclusive), source cut to match. */
MSeg clipTarget(const MSeg &m, int64_t lo, int64_t hi) {
    MSeg r = m;
    r.tLo = lo;
    r.tHi = hi;
    if (m.parallel()) {
        r.sLo = m.sLo + (lo - m.tLo);
        r.sHi = m.sHi - (m.tHi - hi);
    } else {
        r.sLo = m.sLo + (m.tHi - hi);
        r.sHi = m.sHi - (lo - m.tLo);
    }
    return r;
}

/* insertAndBreakOverlaps over the whole set: cut each segment at every boundary of
 * any segment overlapping it on the target, so that overlapping segments become
 * identical on the target; then keep one of each (source, target) pair, ordered. */
vector<MSeg> breakOverlaps(vector<MSeg> segs) {
    vector<MSeg> out;
    sort(segs.begin(), segs.end(), [](const MSeg &a, const MSeg &b) { return a.tLo < b.tLo; });
    size_t i = 0;
    while (i < segs.size()) {
        /* a cluster of segments connected by target overlap */
        size_t j = i + 1;
        int64_t hi = segs[i].tHi;
        while (j < segs.size() && segs[j].tLo <= hi) hi = max(hi, segs[j++].tHi);
        if (j == i + 1) {
            out.push_back(segs[i]);
        } else {
            vector<int64_t> cuts;     /* first positions of pieces */
            for (size_t k = i; k < j; k++) {
                cuts.push_back(segs[k].tLo);
                cuts.push_back(segs[k].tHi + 1);
            }
            sort(cuts.begin(), cuts.end());
            cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());
            for (size_t k = i; k < j; k++) {
                const MSeg &m = segs[k];
                auto c = upper_bound(cuts.begin(), cuts.end(), m.tLo);
                int64_t from = m.tLo;
                for (; c != cuts.end() && *c <= m.tHi; ++c) {
                    out.push_back(clipTarget(m, from, *c - 1));
                    from = *c;
                }
                out.push_back(clipTarget(m, from, m.tHi));
            }
        }
        i = j;
    }
    sort(out.begin(), out.end(), msLess);
    out.erase(unique(out.begin(), out.end(),
                     [](const MSeg &a, const MSeg &b) { return !msLess(a, b) && !msLess(b, a); }),
              out.end());
    return out;
}

/* MappedSegment::canMergeRightWith(next, cutSet = queryCutPoints, sourceCutSet = {}). */
bool canMergeRight(const MSeg &a, const MSeg &b, const set<int64_t> &cutSet) {
    if (a.tRev != b.tRev || a.sRev != b.sRev) return false;
    int64_t qdelta, rdelta, cut;
    if (!a.tRev && !a.sRev) {
        qdelta = b.tStart() - a.tEnd();
        rdelta = b.sStart() - a.sEnd();
        cut = a.tEnd();
    } else if (a.tRev && a.sRev) {
        qdelta = b.tEnd() - a.tStart();
        rdelta = b.sEnd() - a.sStart();
        cut = a.tStart();
    } else if (!a.tRev && a.sRev) {
        qdelta = b.tStart() - a.tEnd();
        rdelta = a.sEnd() - b.sStart();
        cut = a.tEnd();
    } else {
        qdelta = b.tEnd() - a.tStart();
        rdelta = a.sStart() - b.sEnd();
        cut = a.tStart();
    }
    if (qdelta != 1 || rdelta != 1) return false;
    return cutSet.find(cut) == cutSet.end();
}

/* BlockMapper::equalTargetStart */
bool equalTargetStart(const MSeg &a, const MSeg &b) { return a.tLo == b.tLo; }
}

BedLiftover::BedLiftover(const Genome *src, const Genome *tgt, bool dupes, ostream *out)
    : _src(src), _tgt(tgt), _dupes(dupes), _out(out) {}

void BedLiftover::convert(istream *in, int bedType) { scan(in, bedType); }

/* Liftover::visitLine */
void BedLiftover::visitLine() {
    _outBedLines.clear();
    _srcSeq = _src->findSequence(_bedLine._chrName);
    if (_srcSeq < 0) {
        if (_missed.insert(_bedLine._chrName).second)
            cerr << "Unable to find sequence " << _bedLine._chrName << " in genome " << _src->name() << endl;
        return;
    } else if (_bedLine._end > _src->sequenceLength((uint32_t)_srcSeq)) {
        cerr << "Skipping interval with endpoint " << _bedLine._end << "because sequence " << _bedLine._chrName
             << " has length " << _src->sequenceLength((uint32_t)_srcSeq) << endl;
        return;
    } else if (_bedLine._bedType > 9 && _bedLine._blocks.empty()) {
        cerr << "Skipping input line with 0 blocks" << endl;
        return;
    }
    _mappedBlocks.clear();
    if (_bedLine._bedType <= 9) {
        liftBedInterval(_mappedBlocks);
    } else {
        liftBlockIntervals();
    }
    if (!_mappedBlocks.empty() && _bedLine._bedType > 9) assignBlocksToIntervals();
    if (_bedLine._bedType <= 9) _outBedLines = _mappedBlocks;       /* writeBlocksAsIntervals */
    cleanResults();
    _outBedLines.sort(hal::BedLineSrcLess());
    for (hal::BedLine &b : _outBedLines) b.write(*_out);
}

/* BlockLiftover::liftInterval */
void BedLiftover::liftBedInterval(BedList &mappedBedLines) {
    vector<LiftPiece> pieces;
    liftInterval(_src, (uint32_t)_srcSeq, _bedLine._start, _bedLine._end, _tgt, _dupes, pieces, &_cursors);
    if (pieces.empty()) return;
    bool flip = _bedLine._strand == '-';
    int64_t srcOff = _src->sequenceStart((uint32_t)_srcSeq);
    vector<MSeg> segs;
    for (const LiftPiece &p : pieces) {
        int64_t tOff = _tgt->sequenceStart(p.tgtSeq);
        segs.push_back(MSeg{tOff + p.tgtStart, tOff + p.tgtStart + p.length - 1, p.reversed != flip,
                            srcOff + p.srcStart, srcOff + p.srcStart + p.length - 1, flip, p.tgtSeq});
    }
    segs = breakOverlaps(move(segs));

    /* the loop over _mappedSegments with BlockMapper::extractSegment */
    vector<bool> gone(segs.size(), false);
    set<int64_t> queryCutPoints;
    vector<size_t> fragments;
    for (size_t start = 0; start < segs.size(); start++) {
        if (gone[start]) continue;
        fragments.assign(1, start);
        uint32_t startSeq = segs[start].tSeq;
        auto nextAlive = [&](size_t k) {
            while (k < segs.size() && gone[k]) k++;
            return k;
        };
        vector<size_t> v1{start}, v2, toErase;
        size_t next = nextAlive(start + 1);
        while (next < segs.size() && equalTargetStart(segs[v1.back()], segs[next])) {
            v1.push_back(next);
            next = nextAlive(next + 1);
        }
        while (next < segs.size()) {
            while (next < segs.size() && (v2.empty() || equalTargetStart(segs[v2.back()], segs[next])) &&
                   v2.size() < v1.size()) {
                v2.push_back(next);
                next = nextAlive(next + 1);
            }
            bool canMerge = v1.size() == v2.size();
            for (size_t i = 0; i < v1.size() && canMerge; ++i)
                canMerge = segs[v2[i]].tSeq == startSeq && canMergeRight(segs[v1[i]], segs[v2[i]], queryCutPoints);
            if (!canMerge) break;
            fragments.push_back(v2[0]);
            toErase.push_back(v2[0]);
            v1.clear();
            swap(v1, v2);
        }
        if (v1.size() > 1) {
            const MSeg &b = segs[fragments.back()];
            queryCutPoints.insert(max(b.tStart(), b.tEnd()));
        }
        for (size_t e : toErase) gone[e] = true;

        const MSeg &front = segs[fragments.front()], &back = segs[fragments.back()];
        int64_t seqStart = _tgt->sequenceStart(startSeq);
        mappedBedLines.push_back(_bedLine);
        hal::BedLine &out = mappedBedLines.back();
        out._blocks.clear();
        out._chrName = _tgt->sequenceName(startSeq);
        out._start = min(min(front.tStart(), front.tEnd()), min(back.tStart(), back.tEnd())) - seqStart;
        out._end = 1 + max(max(front.tStart(), front.tEnd()), max(back.tStart(), back.tEnd())) - seqStart;
        out._strand = segs[start].tRev ? '-' : '+';
        out._srcStart = min(min(front.sStart(), front.sEnd()), min(back.sStart(), back.sEnd()));
        out._srcStrand = front.sRev ? '-' : '+';
        if (_bedLine._strand == '.') {
            out._strand = '.';
            out._srcStrand = '.';
        }
    }
}

/* Liftover::liftBlockIntervals */
void BedLiftover::liftBlockIntervals() {
    hal::BedLine original = _bedLine;
    std::sort(_bedLine._blocks.begin(), _bedLine._blocks.end());
    for (const hal::BedBlock &blk : _bedLine._blocks) {
        _bedLine._start = blk._start + original._start;
        _bedLine._end = _bedLine._start + blk._length;
        if (_bedLine._end > _bedLine._start) liftBedInterval(_mappedBlocks);
    }
    _bedLine._start = original._start;
    _bedLine._end = original._end;
}

/* Liftover::assignBlocksToIntervals, without PSL */
void BedLiftover::assignBlocksToIntervals() {
    _mappedBlocks.sort(hal::BedLineSrcLess());
    for (auto blockIt = _mappedBlocks.begin(); blockIt != _mappedBlocks.end(); ++blockIt) {
        if (_outBedLines.empty() || !compatible(_outBedLines.back(), *blockIt)) _outBedLines.push_back(*blockIt);
        hal::BedLine &tgtBed = _outBedLines.back();
        tgtBed._start = min(tgtBed._start, blockIt->_start);
        tgtBed._end = max(tgtBed._end, blockIt->_end);
        hal::BedBlock block;
        block._start = blockIt->_start;
        block._length = blockIt->_end - blockIt->_start;
        tgtBed._blocks.push_back(block);
    }
    for (hal::BedLine &b : _outBedLines)
        for (hal::BedBlock &blk : b._blocks) blk._start -= b._start;
    if (!_outBedLines.empty()) flipBlocks(_outBedLines);
}

/* Liftover::compatible */
bool BedLiftover::compatible(const hal::BedLine &tgtBed, const hal::BedLine &newBlock) const {
    if (tgtBed._strand != newBlock._strand) return false;
    if (tgtBed._srcStart == newBlock._srcStart) return false;
    const hal::BedBlock &tgtBlock = tgtBed._blocks.back();
    int64_t delta = tgtBed._strand != _bedLine._strand ? tgtBlock._start - newBlock._end
                                                       : newBlock._start - (tgtBlock._start + tgtBlock._length);
    if (delta < 0) return false;
    return tgtBed._chrName == newBlock._chrName;
}

/* Liftover::flipBlocks, without PSL */
void BedLiftover::flipBlocks(BedList &bedList) const {
    for (hal::BedLine &b : bedList) {
        if (b._blocks.size() > 1) {
            int64_t delta = b._blocks[1]._start - (b._blocks[0]._start + b._blocks[0]._length);
            if (delta < 0) std::reverse(b._blocks.begin(), b._blocks.end());
        }
    }
}

/* Liftover::cleanResults, without PSL */
void BedLiftover::cleanResults() {
    if (_bedLine._bedType <= 6) return;
    for (auto i = _outBedLines.begin(); i != _outBedLines.end();) {
        auto j = next(i);
        if (_bedLine._thickStart != 0 || _bedLine._thickEnd != 0) {
            i->_thickStart = i->_start;
            i->_thickEnd = i->_end;
        }
        if (_bedLine._bedType > 9 && i->_blocks.empty()) _outBedLines.erase(i);
        i = j;
    }
}
}
