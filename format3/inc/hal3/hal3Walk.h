/*
 * HAL format 3 API: every genome's pieces for a reference interval at once.
 *
 * A ReferenceWalk lifts intervals of a reference genome to every printed genome in one
 * pass over the tree: hal maf's columns.  By default a reference base reaches every
 * descendant of the highest ancestral base it goes up to, paralogs included (hal2maf's
 * column); onlyOrthologs lifts to each genome through its common ancestor with the
 * reference instead (halLiftover's).  One ReferenceWalk per thread.
 */
#ifndef _HAL3_API_WALK_H
#define _HAL3_API_WALK_H

#include "hal3Lift.h"

#include <set>
#include <utility>

namespace hal3 {
namespace api {

struct WalkOptions {
    bool ancestors = true;              /* false: no ancestral genome printed (leaf reference) */
    bool dupes = true;                  /* false: canonical copies only */
    bool onlyOrthologs = false;
    std::vector<std::string> targets;   /* only these genomes printed; empty: all */
    std::string rootGenome;             /* only this subtree printed */
    std::vector<std::string> exclude;   /* never printed; the walk still goes through them */
    bool refSingleCopy = false;         /* the reference's other copies not printed */
};

class ReferenceWalk {
  public:
    /* Throws Error(Argument) on unknown genomes or bad combinations. */
    ReferenceWalk(std::shared_ptr<const Alignment> alignment, const Genome *ref, const WalkOptions &opts = WalkOptions());
    ~ReferenceWalk();
    ReferenceWalk(const ReferenceWalk &) = delete;
    ReferenceWalk &operator=(const ReferenceWalk &) = delete;

    /* intervals (sorted, disjoint) of reference sequence refSeq, lifted to every printed
     * genome: pieces[g], source the reference.  The reference's own identity pieces are
     * left out; its other copies are not. */
    void lift(uint32_t refSeq, const std::vector<std::pair<int64_t, int64_t>> &intervals,
              std::map<const Genome *, std::vector<LiftPiece>> &pieces);
    const std::set<const Genome *> &printed() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};
}
}

#endif
