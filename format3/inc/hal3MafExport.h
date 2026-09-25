/*
 * HAL format 3: MAF and TAF export, row by row, built on liftover.
 *
 * For each chunk of a reference sequence, the tree walk (hal3Walk.h) lifts it to every
 * genome in scope at once: hal2maf's columns.
 *
 * Each printed genome's pieces are chained into rows: a piece continues a row when it
 * follows it directly in that genome (same sequence and strand, no bases between) and
 * later on the reference; bases skipped on the reference are deletions, gaps in the
 * row.  A genome that jumps starts a new row; the block is not cut for it.  Blocks
 * are cut where a large share of rows break (breakFrac), where continuing would pad
 * the block with too many gap cells (gapBudget), at the caps, and at chunk and
 * sequence boundaries.
 *
 * Insertions (joinGap, as taffy norm's joining, but done while blocks form): a row
 * also continues across up to joinGap bases of its genome that are in no column with
 * a reference base -- checked exactly, by walking the bases up to their top ancestral
 * bases ("top keys") and down towards the reference.  Those bases get columns without
 * a reference base, placed at the reference boundary where the row's gap is; bases of
 * different rows share such a column when they share a top key, i.e. descend from the
 * same ancestral base, which is exactly HAL's homology (no re-alignment).
 *
 * Rows: the reference row first, then genomes by their distance from the reference in
 * the tree's preorder (cactus-hal2maf's order), and within a genome by sequence name,
 * strand (+ first) and position (taffy sort's).
 *
 * Output, in reference order:
 *   MAF: hal2maf's header, 'a' lines, tab-separated 's' lines named genome.sequence.
 *   TAF: the same blocks as taffy writes them (taffy view -m gives the MAF back): one
 *        line per column, run-length encoded by default, the rows' coordinates at each
 *        block's first column.  A row continues from one block to the next exactly
 *        when it is the same row (no guessing, as taffy must from MAF); a row that
 *        ends and a later one on the same sequence and strand also continue, with the
 *        bases between skipped ("g").  All rows are given coordinates where the
 *        reference's are, as taffy does, so taffy index and view -r work.
 * Either plain or BGZF (blocked gzip, readable by gzip, indexable by taffy index).
 *
 * Chunks and threads: every interval is cut into chunks of opts.chunk bases, and no
 * block, row or insertion crosses a chunk edge, so each chunk is exported on its own:
 * by any of the worker threads, each with its own cursors and DNA handles for the
 * length of the chunk (which bounds memory), and written out in order.  The output is
 * the same, byte for byte, for any number of threads.  No HAL API; DNA comes from the
 * .2bit files.
 */
#ifndef _HAL3MAFEXPORT_H
#define _HAL3MAFEXPORT_H

#include "hal3Alignment.h"
#include "hal3Walk.h"

#include <cstdint>
#include <string>
#include <vector>

namespace hal3 {

struct MafOptions {
    bool ancestors = true;              /* false: --noAncestors */
    bool dupes = true;                  /* false: --noDupes */
    bool onlyOrthologs = false;
    std::vector<std::string> targets;   /* --targetGenomes; empty: all */
    std::string rootGenome;             /* --rootGenome */
    std::vector<std::string> exclude;   /* --excludeGenomes: no rows for these */
    bool refSingleCopy = false;         /* no rows for the reference's other copies */
    double breakFrac = 0.5;
    double gapBudget = 1.0;
    int64_t maxBlockCols = 10000;
    int64_t maxBlockCells = 5000000;
    int64_t chunk = 2000000;
    int64_t joinGap = 30;               /* longest insertion a row continues across; 0: none */
    bool noDna = false;                 /* write N instead of bases */
    bool onlySequenceNames = false;
};

struct MafFormat {
    bool taf = false;
    bool bgzf = false;
    bool runLength = true;              /* TAF: run-length encode each column's bases */
};

struct MafStats {
    uint64_t blocks = 0, rows = 0, columns = 0, alignedBases = 0;
    uint64_t insertionColumns = 0, insertedBases = 0, joinsRejected = 0;
};

/* The blocks of these reference intervals (blocks never cross an interval's ends), in
 * order, to path ("stdout" or a file), exported by `threads` threads.  With one thread
 * each chunk is written as it is made; with more, a chunk's output waits for the chunks
 * before it in memory, or past 32 MB in an unlinked temporary file beside the output
 * (in $TMPDIR, else /tmp, for stdout or when that directory cannot be written).  Throws
 * std::runtime_error. */
MafStats exportMaf(const Alignment &alignment, const Genome *ref, const MafOptions &opts, const MafFormat &format,
                   const std::vector<RefInterval> &intervals, const std::string &path, int threads);
}

#endif
