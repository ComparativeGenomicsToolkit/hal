/*
 * HAL format 3 API: liftover.
 *
 * A Lifter maps intervals between genomes through the tree -- up each branch to the
 * lowest common ancestor, then down -- with every copy of a duplicated base, or with
 * dupes off only the canonical copy (the one with the smallest child coordinate,
 * HAL's canonical paralog in Cactus alignments).  It keeps the chunks it decodes, so
 * nearby queries are cheap: one Lifter per thread.
 */
#ifndef _HAL3_API_LIFT_H
#define _HAL3_API_LIFT_H

#include "hal3Reader.h"

namespace hal3 {
namespace api {

/* [srcStart, srcStart+length) of the source sequence pairs with [tgtStart, ...) of
 * target sequence tgtSeq; reversed pairs srcStart+i with tgtStart+length-1-i. */
struct LiftPiece {
    int64_t srcStart;
    uint32_t tgtSeq;
    int64_t tgtStart;
    int64_t length;
    bool reversed;
};

class Lifter {
  public:
    explicit Lifter(std::shared_ptr<const Alignment> alignment);
    ~Lifter();
    Lifter(const Lifter &) = delete;
    Lifter &operator=(const Lifter &) = delete;

    /* [start, end) of sequence srcSeq of src, lifted to tgt: out replaced with the
     * aligned pieces */
    void lift(const Genome *src, uint32_t srcSeq, int64_t start, int64_t end, const Genome *tgt, bool dupes,
              std::vector<LiftPiece> &out);
    /* One branch at a time, for callers that walk the tree themselves.  A piece relates
     * a source interval to the genome reached so far.  up: from g to its parent, the
     * parts with no parent to *unmapped (canonicalOnly: also the copies that are not
     * canonical).  down: from child's parent to child. */
    void up(const Genome *g, const std::vector<LiftPiece> &in, std::vector<LiftPiece> &out,
            std::vector<LiftPiece> *unmapped = nullptr, bool canonicalOnly = false);
    void down(const Genome *child, bool dupes, const std::vector<LiftPiece> &in, std::vector<LiftPiece> &out);

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};
}
}

#endif
