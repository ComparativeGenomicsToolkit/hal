/*
 * HAL format 3 API: opening an alignment, its tree and its genomes' sequences.
 *
 * An Alignment is immutable and may be shared by threads; open() returns it as a
 * shared_ptr, which the per-thread readers (Lifter, ReferenceWalk) keep alive.  The
 * tree, genome lengths and sequence counts come from the manifest, so listing an
 * alignment opens no genome file; a genome's sequence table is read the first time it
 * is asked for, and kept.
 *
 * Coordinates: a position is (sequence index, position in the sequence), 0-based and
 * half-open.  Genome coordinates -- positions in the genome's sequences laid end to
 * end, as HAL has them -- are available through sequenceStart and sequenceAt.
 */
#ifndef _HAL3_API_READER_H
#define _HAL3_API_READER_H

#include "hal3Error.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace hal3 {
namespace api {

class Genome;

struct OpenOptions {
    /* A local directory's members share at most this many open files (0: the soft limit
     * on open files / 4, at most 256; HAL3_MAX_OPEN_FILES overrides it).  A pack is one
     * open file. */
    size_t maxOpenFiles = 0;
};

class Alignment {
  public:
    /* A pack file, a directory or a manifest; a local path or an http(s) URL. */
    static std::shared_ptr<const Alignment> open(const std::string &path, const OpenOptions &opts = OpenOptions());
    ~Alignment();
    Alignment(const Alignment &) = delete;
    Alignment &operator=(const Alignment &) = delete;

    const std::string &path() const;
    bool isPack() const;
    const Genome *root() const;
    const Genome *genome(const std::string &name) const;   /* nullptr if absent */
    /* every genome, the root first, then the manifest's order (a preorder) */
    const std::vector<const Genome *> &genomes() const;
    const std::string &newick() const;
    const std::map<std::string, std::string> &meta() const;

    struct Impl;   /* the library's own: opaque here */
    const Impl &impl() const { return *_impl; }

  private:
    Alignment();
    std::unique_ptr<Impl> _impl;
};

class Genome {
  public:
    const std::string &name() const;
    const Genome *parent() const;                            /* nullptr for the root */
    const std::vector<const Genome *> &children() const;
    int depth() const;                                       /* the root is 0 */
    double branchLength() const;                             /* to the parent; -1 for the root */
    int64_t length() const;                                  /* bases */
    int64_t numSequences() const;
    bool hasDna() const;
    std::vector<std::string> memberRoles() const;            /* extra members, e.g. "lik" */

    /* the sequence table (read on first use) */
    const std::string &sequenceName(uint32_t seq) const;
    int64_t sequenceLength(uint32_t seq) const;
    int64_t sequenceStart(uint32_t seq) const;               /* genome coordinate of its base 0 */
    int64_t findSequence(const std::string &name) const;     /* -1 if absent */
    uint32_t sequenceAt(int64_t genomePos) const;            /* the sequence holding it */
    /* the genome's metadata (HAL's per-genome metadata, as converted) */
    std::map<std::string, std::string> meta() const;

    struct Impl;   /* the library's own: opaque here */
    const Impl &impl() const { return *_impl; }
    ~Genome();
    Genome(const Genome &) = delete;
    Genome &operator=(const Genome &) = delete;

  private:
    friend class Alignment;
    Genome();
    std::unique_ptr<Impl> _impl;
};

/* The lowest common ancestor of two genomes of one alignment. */
const Genome *lowestCommonAncestor(const Genome *a, const Genome *b);
}
}

#endif
