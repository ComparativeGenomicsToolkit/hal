/*
 * HAL format 3: MAF export, row by row, built on liftover.
 *
 * For each chunk of a reference sequence, one tree walk lifts it to every genome in
 * scope: up the reference's lineage to the scope top, collecting at each genome the
 * bases that go no higher ("tops"), then one pass down from the top that descends from
 * every top into its whole subtree.  That is hal2maf's column: every descendant of the
 * highest ancestral base the reference base reaches, paralogs included.
 * --onlyOrthologs instead descends at each ancestor only into the children off the
 * reference's lineage (and into the reference's own subtree), which is liftover to the
 * pairwise common ancestor; --noDupes goes up only from canonical copies and down only
 * to them.
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
 * Output is MAF (hal2maf's header, 'a' lines, tab-separated 's' lines named
 * genome.sequence), in reference order; the reference row first, then rows by genome
 * name, sequence and position.  No HAL API; DNA comes from the .2bit files.
 */
#ifndef _HAL3MAFEXPORT_H
#define _HAL3MAFEXPORT_H

#include "hal3Alignment.h"
#include "hal3Liftover.h"
#include "twobit64.h"

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace hal3 {

struct MafOptions {
    bool ancestors = true;              /* false: --noAncestors */
    bool dupes = true;                  /* false: --noDupes */
    bool onlyOrthologs = false;
    std::vector<std::string> targets;   /* --targetGenomes; empty: all */
    std::string rootGenome;             /* --rootGenome */
    double breakFrac = 0.5;
    double gapBudget = 1.0;
    int64_t maxBlockCols = 10000;
    int64_t maxBlockCells = 5000000;
    int64_t chunk = 2000000;
    int64_t joinGap = 30;               /* longest insertion a row continues across; 0: none */
    bool noDna = false;                 /* write N instead of bases */
    bool onlySequenceNames = false;
};

/* Where MAF text goes: plain, or BGZF (blocked gzip, readable by gzip and indexable). */
class MafSink {
  public:
    virtual ~MafSink() {}
    virtual void write(const std::string &text) = 0;
    virtual void close() = 0;
};
/* "stdout" or a path.  Throws std::runtime_error. */
std::unique_ptr<MafSink> openMafSink(const std::string &path, bool bgzf);

class MafExporter {
  public:
    MafExporter(const Alignment &alignment, const Genome *ref, const MafOptions &opts, MafSink &out);
    ~MafExporter();
    MafExporter(const MafExporter &) = delete;
    MafExporter &operator=(const MafExporter &) = delete;

    void writeHeader();
    /* [start, end) of one reference sequence; blocks never cross its ends. */
    void exportInterval(uint32_t refSeq, int64_t start, int64_t end);

    struct Stats {
        uint64_t blocks = 0, rows = 0, columns = 0, alignedBases = 0;
        uint64_t insertionColumns = 0, insertedBases = 0, joinsRejected = 0;
    };
    const Stats &stats() const { return _stats; }

  private:
    struct Chain;
    struct Key {   /* a top ancestral base: where a base's homology is decided */
        const Genome *node;
        uint32_t seq;
        int64_t pos;
        bool operator<(const Key &o) const;
        bool operator==(const Key &o) const { return node == o.node && seq == o.seq && pos == o.pos; }
    };
    const Alignment &_aln;
    const Genome *_ref;
    MafOptions _opts;
    MafSink &_out;
    std::set<const Genome *> _scope, _printed;
    std::vector<const Genome *> _refPath;   /* the reference, its parent, ... up to _top */
    const Genome *_top = nullptr;
    CursorCache _cursors;
    std::map<const Genome *, tb_file *> _dna;
    std::map<const Genome *, std::vector<int64_t>> _dnaIndex;
    Stats _stats;

    void chunk(uint32_t refSeq, int64_t start, int64_t end);
    void walk(uint32_t refSeq, int64_t start, int64_t end, std::map<const Genome *, std::vector<LiftPiece>> &rows);
    void descend(const Genome *v, std::vector<LiftPiece> &in, const Genome *skipChild,
                 std::map<const Genome *, std::vector<LiftPiece>> &tops,
                 std::map<const Genome *, std::vector<LiftPiece>> &rows);
    bool insertionKeys(const Genome *g, uint32_t seq, bool reversed, int64_t o1, int64_t o2, std::vector<Key> &keys);
    void fetch(const Genome *g, uint32_t seq, int64_t start, int64_t end, bool reversed, std::string &out);
};
}

#endif
