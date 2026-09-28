/*
 * HAL format 3: liftover.
 *
 * Two layers.
 *
 * liftInterval maps an interval of one genome to another through the tree: up each
 * branch to the lowest common ancestor with an axis-a query on the genome's .hgn (a
 * base has at most one parent), then down each branch with an axis-b query on the
 * child's .hgn (every child copy of a parent base -- paralogy -- or with dupes off only
 * the copy with the smallest child genome coordinate, HAL's canonical paralog).  The
 * result is the aligned base pairs, as pieces.  These are HAL's mapUp / mapDown /
 * mapSelf semantics, on maximal runs rather than segments.
 *
 * Writing BED as halLiftover does is in hal3BedLiftover.h.
 */
#ifndef _HAL3LIFTOVER_H
#define _HAL3LIFTOVER_H

#include "hal3Alignment.h"

#include <list>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace hal3 {

/* [srcStart, srcStart+length) of the source sequence pairs with [tgtStart, ...) of
 * target sequence tgtSeq; reversed pairs srcStart+i with tgtStart+length-1-i. */
struct LiftPiece {
    int64_t srcStart;
    uint32_t tgtSeq;
    int64_t tgtStart;
    int64_t length;
    bool reversed;
};

/* Open liftmap cursors kept between queries.  Each caches the chunks it decodes, so
 * nearby queries -- the pieces of one interval, neighbouring BED lines -- decode each
 * chunk once instead of once per query.  Least recently used cursors are closed past
 * maxCursors.  Not thread-safe: one per thread. */
class CursorCache {
  public:
    explicit CursorCache(size_t bytesPerCursor = (size_t)8 << 20, size_t maxCursors = 64);
    ~CursorCache();
    CursorCache(const CursorCache &) = delete;
    CursorCache &operator=(const CursorCache &) = delete;
    /* The runs of member m of `axis` of g's .hgn overlapping [lo, hi), clipped to it. */
    void query(const Genome *g, int axis, uint32_t m, int64_t lo, int64_t hi, std::vector<lmap_hit> &out);

  private:
    typedef std::tuple<lmap_file *, int, uint32_t> Key;
    struct KeyHash {
        size_t operator()(const Key &k) const {
            return std::hash<const void *>()(std::get<0>(k)) ^ ((size_t)std::get<2>(k) * 2 + (size_t)std::get<1>(k)) * 0x9e3779b97f4a7c15ULL;
        }
    };
    struct Entry {
        lmap_cursor *cursor;
        std::list<Key>::iterator use;   /* its place in _lru */
    };
    size_t _bytes, _max;
    std::list<Key> _lru;                /* most recently used first */
    std::unordered_map<Key, Entry, KeyHash> _open;
};

/* One branch at a time.  A piece relates a source interval to an interval of the
 * genome currently reached (tgtSeq, tgtStart); these carry pieces over one branch.
 *
 * liftUpBranch: from g to its parent.  The parts with no parent position go to
 * *unmapped (still in g's coordinates) when it is given.  canonicalOnly (--noDupes
 * walks): a base goes up only if it is the canonical copy of its parent base (the copy
 * with the smallest child coordinate); other copies stay behind as unmapped.
 *
 * liftDownBranch: from child's parent to child, every copy, or with dupes off the
 * canonical copy of each parent base. */
void liftUpBranch(const Genome *g, const std::vector<LiftPiece> &in, std::vector<LiftPiece> &out,
                  std::vector<LiftPiece> *unmapped, CursorCache &cursors, bool canonicalOnly = false);
void liftDownBranch(const Genome *child, bool dupes, const std::vector<LiftPiece> &in, std::vector<LiftPiece> &out,
                    CursorCache &cursors);

/* Map [start, end) of sequence srcSeq of src to tgt.  Throws std::runtime_error. */
void liftInterval(const Genome *src, uint32_t srcSeq, int64_t start, int64_t end, const Genome *tgt, bool dupes,
                  std::vector<LiftPiece> &out, CursorCache *cursors = nullptr);

}

#endif
