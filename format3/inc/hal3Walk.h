/*
 * HAL format 3: the tree walk behind reference-anchored output (hal maf, hal coverage).
 *
 * Lifts intervals of a reference genome to every genome in scope at once.
 *
 * Default: up the reference's lineage to the scope top, collecting at each genome the
 * bases that go no higher ("tops"), then one pass down from the top that descends from
 * every top into its whole subtree.  That is hal2maf's column: every descendant of the
 * highest ancestral base a reference base reaches, paralogs included.
 *
 * onlyOrthologs: at each ancestor on the reference's lineage, descend only into the
 * children off the lineage (and into the reference's own subtree): liftover to each
 * genome's common ancestor with the reference, as halLiftover and halCoverage.
 *
 * dupes off (--noDupes): up only from canonical copies, down only to them.
 *
 * Scope: the spanning tree of the printed genomes -- the targets, or the subtree of
 * rootGenome, or every genome -- and the reference.  Only printed genomes get pieces;
 * ancestors are left out with ancestors = false.
 */
#ifndef _HAL3WALK_H
#define _HAL3WALK_H

#include "hal3Alignment.h"
#include "hal3Liftover.h"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace hal3 {

struct WalkOptions {
    bool ancestors = true;
    bool dupes = true;
    bool onlyOrthologs = false;
    std::vector<std::string> targets;   /* empty: every genome */
    std::string rootGenome;
};

class TreeWalk {
  public:
    /* Throws std::runtime_error on unknown genomes or bad combinations. */
    TreeWalk(const Alignment &alignment, const Genome *ref, const WalkOptions &opts);

    /* Lift intervals of the reference sequence refSeq (sorted, disjoint) to every
     * printed genome.  pieces[g]: source = the reference, target = g.  The reference's
     * own identity pieces are left out; its paralogous copies are not. */
    void lift(uint32_t refSeq, const std::vector<std::pair<int64_t, int64_t>> &intervals,
              std::map<const Genome *, std::vector<LiftPiece>> &pieces);

    const Genome *ref() const { return _ref; }
    const Genome *top() const { return _top; }
    /* the reference, its parent, ... up to top() */
    const std::vector<const Genome *> &refPath() const { return _refPath; }
    const std::set<const Genome *> &printed() const { return _printed; }
    const WalkOptions &options() const { return _opts; }
    CursorCache &cursors() { return _cursors; }

  private:
    const Alignment &_aln;
    const Genome *_ref;
    WalkOptions _opts;
    std::set<const Genome *> _scope, _printed;
    std::vector<const Genome *> _refPath;
    const Genome *_top = nullptr;
    CursorCache _cursors;

    void descend(const Genome *v, std::vector<LiftPiece> &in, std::map<const Genome *, std::vector<LiftPiece>> &tops,
                 std::map<const Genome *, std::vector<LiftPiece>> &pieces);
};

/* Merge pieces contiguous in both source and target. */
void coalescePieces(std::vector<LiftPiece> &v);
}

#endif
