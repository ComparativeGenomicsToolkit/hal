/*
 * HAL format 3: universal columns as trees of bases, for scoring them where they are made
 * (phyloU's native input) instead of through a MAF.
 *
 * The columns are hal maf --universal's with duplicates kept -- the same index, planner
 * and axis (hal3Universal.h), so hal wig projects scores on them -- and every base of
 * every genome in scope is in exactly one, ancestors' included.  Each base also has its
 * parent: the base of the parent genome it aligns to, which is always in the same column.
 * A column is so a tree of bases -- the species tree with a genome repeated once per copy
 * (a copy tree) -- rooted at its one base with no parent, the column's top.
 *
 * There are no rows and no blocks: a chunk is its columns, each a list of entries in
 * genome preorder (a parent before its children).  An entry's piece is the stretch of its
 * genome it lies in, carried down from one piece of the parent genome: two entries of a
 * genome under one parent entry are copies of it, and the same two pieces under one
 * parent piece over several columns are the same duplication.
 */
#ifndef _HAL3COLUMNS_H
#define _HAL3COLUMNS_H

#include "hal3Alignment.h"
#include "hal3MafExport.h"
#include "hal3Universal.h"
#include "hal3Walk.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace hal3 {

struct ColumnOptions {
    std::string rootGenome;              /* the backbone; "" for the root */
    std::vector<std::string> exclude;    /* leaves left out */
    bool leafInsertions = false;         /* a leaf's bases with no parent: columns of their own */
    bool ancestorBases = true;           /* false: ancestors' bases are N, their DNA not read */
    int64_t chunk = (int64_t)1 << 20;    /* columns a chunk */
};

/* Columns [u, u+length) of a chunk hold bases [start, start+length) of sequence seq of
 * genome `genome`, in increasing order, or in decreasing order when reversed (and then
 * complemented). */
struct ColumnPiece {
    uint32_t genome, seq;
    int64_t start, length, u;
    bool reversed;
};

struct ColumnChunk {
    const std::string *axisName = nullptr;   /* the axis sequence the chunk lies on */
    int64_t axisStart = 0;                   /* where on it the chunk's column 0 is */
    int64_t cols = 0;
    std::vector<int64_t> colStart;   /* cols + 1: column u's entries are [colStart[u], colStart[u+1]) */
    std::vector<uint32_t> genome;    /* per entry: an index into UniversalColumns' genomes */
    std::vector<int32_t> parent;     /* per entry: its parent's entry, counted from colStart[u]; -1: a top */
    std::vector<uint32_t> piece;     /* per entry: an index into pieces */
    std::string base;                /* per entry, in the axis's orientation, soft-masking kept */
    std::vector<ColumnPiece> pieces;
};

class UniversalColumns {
  public:
    /* The columns of these backbone ranges (sorted, disjoint), or of the whole alignment
     * when ranges is empty.  workers: how many threads will call chunk() at once (and
     * build the index).  Throws std::runtime_error. */
    UniversalColumns(const Alignment &alignment, const std::vector<RefInterval> &ranges, const ColumnOptions &opts,
                     int workers);
    ~UniversalColumns();
    UniversalColumns(const UniversalColumns &) = delete;
    UniversalColumns &operator=(const UniversalColumns &) = delete;

    size_t chunks() const { return _plan->chunks(); }
    int64_t columns() const { return _plan->columns(); }
    /* Chunk k by worker w (0 <= w < workers): different workers may run at once.
     * Throws std::runtime_error. */
    void chunk(size_t k, int w, ColumnChunk &out);

    /* the genomes in scope, in preorder from the backbone: entries' genome indices */
    uint32_t genomes() const { return (uint32_t)_genomes.size(); }
    const Genome *genome(uint32_t i) const { return _genomes[i]; }
    int32_t parentOf(uint32_t i) const { return _parentOf[i]; }   /* -1: the backbone */
    const Genome *backbone() const { return _genomes[0]; }

  private:
    struct Worker;
    const Alignment &_aln;
    ColumnOptions _opts;
    WalkOptions _walkOpts;
    std::vector<const Genome *> _genomes;
    std::vector<int32_t> _parentOf;
    std::map<const Genome *, uint32_t> _indexOf;
    std::unique_ptr<UniversalIndex> _index;
    std::unique_ptr<UniversalPlanner> _plan;
    std::vector<AxisRange> _axes;
    std::vector<std::unique_ptr<Worker>> _workers;
};
}

#endif
