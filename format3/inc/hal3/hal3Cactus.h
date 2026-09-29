/*
 * HAL format 3 API: a Cactus subproblem, as cactus_consolidated describes it, to a
 * fragment (hal3Writer.h) -- the replacement for halAppendCactusSubtree.
 *
 * A subproblem aligns an ancestor with its children (and outgroups, which are not kept).
 * Cactus gives it in the .c2h form: for each sequence, in order, a line naming its
 * genome ("event"), its name and whether it is the ancestor's (a "bottom" sequence),
 * then its segments, which tile it:
 *
 *   ancestor:  a <name> <start> <length>                 (named, for children to refer to)
 *   child:     a <start> <length> <ancestor segment> <strand>   (aligned; strand 1:
 *                                                        the same, 0: reversed)
 *              a <start> <length>                        (unaligned)
 *
 * The ancestor's sequences come first.  Each sequence's DNA is given with it.
 *
 * The fragment it becomes: the ancestor as its root (sequences and DNA, no branch), and
 * each child with its branch to the ancestor -- a run per aligned segment -- and its DNA.
 * Merged with its parent's fragment (where the ancestor is a child) and its children's
 * (where they are roots), it makes the alignment; so an ancestor's sequences must be in
 * the same order in both .c2h files it is in, which is Cactus's.
 */
#ifndef _HAL3_API_CACTUS_H
#define _HAL3_API_CACTUS_H

#include "hal3Writer.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace hal3 {
namespace api {

class C2hImporter {
  public:
    /* The subproblem's tree -- the ancestor and its (ingroup) children, with branch
     * lengths, as Cactus writes it for halAppendCactusSubtree -- and the fragment's
     * directory, which must not exist or hold an alignment.  Genomes the tree does not
     * name (outgroups) are skipped. */
    C2hImporter(const std::string &fragmentDir, const std::string &newick, const GenomeWriterOptions &opts = GenomeWriterOptions());
    ~C2hImporter();
    C2hImporter(const C2hImporter &) = delete;
    C2hImporter &operator=(const C2hImporter &) = delete;

    /* First pass: every sequence, in .c2h order, with its length (its segments' total). */
    void declare(const std::string &genome, const std::string &name, bool ancestor, int64_t length);
    /* Second pass, in the same order: each sequence's segment lines (the .c2h text after
     * its "s" line, up to the next one) and its DNA.  Throws Error(Format) on a malformed
     * line, a segment outside its sequence, or a reference to an unknown ancestor
     * segment. */
    void sequence(const std::string &genome, const std::string &name, const char *segments, size_t segmentsLength,
                  const char *dna, int64_t dnaLength);
    /* Close the genomes and write the fragment's manifest (then validated). */
    void finish();

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

/* The same from files: the .c2h, its FASTA (a record per .c2h sequence, in the same
 * order, named as the sequence), and the subproblem's tree (a Newick file or string). */
void importC2h(const std::string &c2hPath, const std::string &fastaPath, const std::string &newick, const std::string &fragmentDir,
               const GenomeWriterOptions &opts = GenomeWriterOptions());
}
}

#endif
