/*
 * HAL format 3 API: the runs of a branch -- HAL's segments, replaced.
 *
 * The branch above a genome (the child) is a set of runs: maximal gapless, colinear
 * stretches of the child aligned to its parent.  A child base is in at most one run; a
 * parent base may be in several (its copies in the child: paralogy).  With dupes off,
 * a parent base's canonical copy is the one with the smallest child coordinate (HAL's
 * canonical paralog in Cactus alignments).  A BranchReader keeps the chunks it decodes,
 * so nearby queries are cheap: one per thread.
 */
#ifndef _HAL3_API_BRANCH_H
#define _HAL3_API_BRANCH_H

#include "hal3Reader.h"

#include <functional>

namespace hal3 {
namespace api {

/* [childStart, +length) of child sequence childSeq pairs with [parentStart, +length) of
 * parent sequence parentSeq; reversed pairs childStart+i with parentStart+length-1-i. */
struct Run {
    uint32_t childSeq, parentSeq;
    int64_t childStart, parentStart, length;
    bool reversed;
};

class BranchReader {
  public:
    explicit BranchReader(std::shared_ptr<const Alignment> alignment);
    ~BranchReader();
    BranchReader(const BranchReader &) = delete;
    BranchReader &operator=(const BranchReader &) = delete;

    /* the runs of the branch above child overlapping [lo, hi) of its sequence childSeq,
     * in child order, clipped to it (clip false: whole runs) */
    void byChild(const Genome *child, uint32_t childSeq, int64_t lo, int64_t hi, std::vector<Run> &out, bool clip = true);
    /* the runs overlapping [lo, hi) of sequence parentSeq of child's parent, in parent
     * order: every copy, clipped unless clip is false; or with dupes off only the
     * canonical copy of each parent base, cut to where it is canonical */
    void byParent(const Genome *child, uint32_t parentSeq, int64_t lo, int64_t hi, bool dupes, std::vector<Run> &out,
                  bool clip = true);
    /* every run of the branch above child, in child order, streamed (nothing kept):
     * whole-branch scans */
    void scan(const Genome *child, const std::function<void(const Run &)> &visit);
    /* the parent base of base pos of child sequence seq; false if it has none */
    bool parentOf(const Genome *child, uint32_t seq, int64_t pos, uint32_t &parentSeq, int64_t &parentPos, bool &reversed);

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};
}
}

#endif
