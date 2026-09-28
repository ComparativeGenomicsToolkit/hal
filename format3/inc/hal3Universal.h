/*
 * HAL format 3: the universal column axis, for hal maf --universal.
 *
 * Scope: the backbone genome B (the root, or --rootGenome) and its subtree, ranked in
 * the tree's preorder (the Newick's), B first.
 *
 * Columns.  A top base is a base with no parent: every base of B, and in any other
 * genome a base no run of its branch covers (with dupes off, also a base that is not
 * the canonical copy of its parent base).  A top interval J = (G, q, [s, e)) is a
 * maximal run of top bases of sequence q of a genome G below B.  A column is the
 * down-closure of one top base: it and every copy of it below, as the tree walk finds
 * them.  Every base of B is a column, and so is every base of every covered J -- an
 * internal genome's J always, a leaf's only with leafInsertions and never an excluded
 * leaf's.  A base has one path up, so it is in exactly one column when its top base is
 * covered.
 *
 * Order ("flank1").  J goes right after its left flank L = s-1 in G's reading
 * direction, else right before its right flank R = e, else (a whole sequence with no
 * aligned base) after the backbone, as a range of its own.  The flank's top base K, a
 * base of a strict ancestor A of G, and the flank's strand relative to A pick a slot of
 * A next to K; the intervals of one slot go newest nearest their anchor base.  The axis
 * is the depth-first expansion of the top-level ranges: each run's bases in order, with
 * each slot's intervals expanded in place, a run reversed when its interval's strand
 * relative to A differs from the orientation of the run holding K.
 *
 * UniversalIndex holds every covered J by slot, built in parallel, one flank lift per J.
 * UniversalPlanner walks it into the axis, cut into chunks, any of which can be walked
 * again on its own from a snapshot.
 */
#ifndef _HAL3UNIVERSAL_H
#define _HAL3UNIVERSAL_H

#include "hal3Alignment.h"
#include "hal3Liftover.h"

#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace hal3 {

/* Positions of genome g are coded code(g, q, k) = sequenceStart(q) + q + k, for k in
 * [0, length(q)]: base k of sequence q and the slot just before it share a code, and
 * successive sequences' codes neither overlap nor touch. */

/* One covered top interval J, placed. */
struct Placement {
    int64_t slot;    /* code of its slot in the anchor genome A; -1: flankless */
    int64_t jcode;   /* code of J's first base in G */
    uint64_t bits;   /* length << 24 | rank of G << 2 | REV << 1 | RIGHT */
    int64_t len() const { return (int64_t)(bits >> 24); }
    int genome() const { return (int)((bits >> 2) & 0x3fffff); }
    bool rev() const { return (bits & 2) != 0; }     /* J reads against A's forward strand */
    bool right() const { return (bits & 1) != 0; }   /* the anchor base is right of the slot */
    int64_t anchor() const { return right() ? slot : slot - 1; }   /* code of the anchor base */
    static uint64_t pack(int64_t len, int genome, bool rev, bool right) {
        return (uint64_t)len << 24 | (uint64_t)genome << 2 | (uint64_t)rev << 1 | (uint64_t)right;
    }
};

/* The order of placements in a slot, in A's forward direction: those anchored left of
 * it first, by (rank, jcode) descending; then those anchored right of it, ascending --
 * newest nearest its anchor base. */
bool placementLess(const Placement &x, const Placement &y);

/* Count, sum and sum of squares of a set of genome coordinates (sequenceStart + pos),
 * mod 2^64: three numbers that catch a base missed, doubled or misplaced. */
struct Moments {
    uint64_t n = 0, s1 = 0, s2 = 0;
    void add(int64_t a, int64_t b);   /* every coordinate in [a, b) */
    void add(const Moments &m) {
        n += m.n;
        s1 += m.s1;
        s2 += m.s2;
    }
    void subtract(const Moments &m) {
        n -= m.n;
        s1 -= m.s1;
        s2 -= m.s2;
    }
    bool operator==(const Moments &m) const { return n == m.n && s1 == m.s1 && s2 == m.s2; }
};

class UniversalIndex {
  public:
    /* The top intervals of every genome below backbone (covered: placed; with allTops,
     * the others are counted too), on `threads` threads.  Throws std::runtime_error. */
    UniversalIndex(const Alignment &alignment, const Genome *backbone, bool leafInsertions, bool dupes,
                   const std::set<const Genome *> &excluded, bool allTops, int threads);

    int genomes() const { return (int)_genomes.size(); }
    const Genome *genome(int rank) const { return _genomes[(size_t)rank]; }
    int rank(const Genome *g) const;   /* -1 outside the scope */
    bool covered(int rank) const { return _covered[(size_t)rank] != 0; }
    int64_t code(int rank, uint32_t q, int64_t k) const { return _codeStart[(size_t)rank][q] + k; }
    void decode(int rank, int64_t code, uint32_t &q, int64_t &k) const;

    /* the placements anchored in a genome, by slot then placementLess */
    const std::vector<Placement> &anchoredIn(int rank) const { return _anchored[(size_t)rank]; }
    /* whole sequences with no aligned base, by (rank, jcode) */
    const std::vector<Placement> &flankless() const { return _flankless; }
    /* a genome's top intervals and their bases (the backbone's: its sequences) */
    uint64_t topIntervals(int rank) const { return _topIntervals[(size_t)rank]; }
    const Moments &topBases(int rank) const { return _topBases[(size_t)rank]; }
    const Moments &allBases(int rank) const { return _allBases[(size_t)rank]; }

    uint64_t placements() const { return _placements; }   /* covered intervals, flankless ones included */
    int64_t columns() const { return _columns; }          /* the backbone's bases plus the covered intervals' */

  private:
    std::vector<const Genome *> _genomes;
    std::map<const Genome *, int> _rank;
    std::vector<char> _covered;
    std::vector<std::vector<int64_t>> _codeStart;   /* per genome, per sequence */
    std::vector<std::vector<Placement>> _anchored;
    std::vector<Placement> _flankless;
    std::vector<uint64_t> _topIntervals;
    std::vector<Moments> _topBases, _allBases;
    uint64_t _placements = 0;
    int64_t _columns = 0;
};

/* A top-level range: bases [start, end) of sequence seq of the genome ranked g. */
struct TopRange {
    int g;
    uint32_t seq;
    int64_t start, end;
};

/* Columns [u, u+len) of a chunk (u counts from its first column): bases [pos, pos+len)
 * of sequence seq of the genome ranked g, in increasing order, or decreasing when rev.
 * depth 0: the top-level range's own bases; d: an interval nested d deep. */
struct Segment {
    int64_t u;
    int g;
    uint32_t seq;
    int64_t pos, len;
    bool rev;
    int depth;
};

class UniversalPlanner {
  public:
    /* backbone: ranges of backbone sequences, in output order, disjoint.  whole: they are
     * every backbone base, so the flankless intervals follow them as ranges of their own
     * and run() checks that every placement and column is visited once.  Each range is
     * cut into chunks of `chunk` columns from its own first column. */
    UniversalPlanner(const UniversalIndex &index, const std::vector<TopRange> &backbone, bool whole, int64_t chunk);

    /* Walk the whole axis once, keeping a snapshot at every chunk.  Throws
     * std::runtime_error if the totals do not add up (whole). */
    void run();

    struct Chunk {
        size_t range;     /* index into ranges() */
        int64_t u0;       /* its first column */
        int64_t cols;
        int64_t rangePos; /* position, in the range's sequence, of its first own base */
    };
    const std::vector<TopRange> &ranges() const { return _ranges; }
    size_t chunks() const { return _chunks.size(); }
    const Chunk &chunk(size_t k) const { return _chunks[k]; }
    /* chunk k's segments, in axis order; from its snapshot, so any thread, any order */
    void segments(size_t k, std::vector<Segment> &out) const;

    int64_t columns() const { return _u; }
    uint64_t visited() const { return _visited; }
    int maxDepth() const { return _maxDepth; }

  private:
    struct Frame {
        int g, depth;
        bool rev;
        int64_t lo, hi;   /* codes of the run's bases, [lo, hi) */
        int64_t cur;      /* forward: the next base; reversed: one past it */
        size_t b, e, it;  /* anchoredIn(g)[b, e): the run's placements; it: the next (reversed: one past it) */
    };
    const UniversalIndex &_idx;
    std::vector<TopRange> _ranges;
    bool _whole;
    int64_t _chunk;
    std::vector<Chunk> _chunks;
    std::vector<std::vector<Frame>> _snaps;
    int64_t _u = 0;
    uint64_t _visited = 0;
    int _maxDepth = 0;

    Frame frame(int g, int64_t lo, int64_t hi, bool rev, int depth) const;
    int64_t step(std::vector<Frame> &st, int64_t want, std::vector<Segment> *out, uint64_t &visited, int &maxDepth) const;
};

/* The top intervals of sequence seq of g, in order, as emit(start, end): the gaps
 * between the runs of its branch, or with canonicalOnly (--noDupes) also its bases
 * that are not the canonical copy of their parent base. */
void forEachTopInterval(const Genome *g, uint32_t seq, bool canonicalOnly, CursorCache &cursors,
                        const std::function<void(int64_t, int64_t)> &emit);
}

#endif
