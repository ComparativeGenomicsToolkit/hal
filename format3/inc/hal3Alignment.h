/*
 * HAL format 3: reading an alignment.
 *
 * An Alignment is the manifest plus the genome files it names, opened as they are
 * first used.  A Genome knows its place in the tree, its sequence table (axis a of
 * its .hgn) and its .hgn, which holds the branch to its parent.  Depends on liftmap
 * only -- no HAL API, no HDF5.
 *
 * Thread safety: an Alignment may be shared by threads; a genome's .hgn is opened
 * once, under a lock, and liftmap files are safe to share.
 */
#ifndef _HAL3ALIGNMENT_H
#define _HAL3ALIGNMENT_H

#include "hal3Manifest.h"
#include "liftmap.h"

#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace hal3 {

class Alignment;

class Genome {
  public:
    const std::string &name() const { return _entry->name; }
    const Genome *parent() const { return _parent; }
    const std::vector<const Genome *> &children() const { return _children; }
    int depth() const { return _depth; }        /* the root is 0 */

    /* The genome file: axis a is this genome's sequences, axis b its parent's.
     * Opened on first use; throws std::runtime_error if it cannot be. */
    lmap_file *hgn() const;

    /* The sequence table, in genome order.  Genome coordinates are sequence
     * coordinates offset by sequenceStart, as in HAL. */
    uint32_t numSequences() const;
    const char *sequenceName(uint32_t i) const;
    int64_t sequenceLength(uint32_t i) const;
    int64_t sequenceStart(uint32_t i) const;
    int64_t findSequence(const std::string &name) const;   /* -1 if absent */

  private:
    friend class Alignment;
    const Alignment *_alignment = nullptr;
    const GenomeEntry *_entry = nullptr;
    const Genome *_parent = nullptr;
    std::vector<const Genome *> _children;
    int _depth = 0;
    mutable std::once_flag _opened;
    mutable lmap_file *_hgn = nullptr;
    mutable std::vector<int64_t> _starts;      /* numSequences + 1 */
    void open() const;
};

class Alignment {
  public:
    /* A format 3 alignment directory (or URL).  Throws std::runtime_error. */
    explicit Alignment(const std::string &dir);
    ~Alignment();
    Alignment(const Alignment &) = delete;
    Alignment &operator=(const Alignment &) = delete;

    const Manifest &manifest() const { return _manifest; }
    const std::string &dir() const { return _dir; }
    const Genome *root() const { return _genomes.empty() ? nullptr : _genomes[0].get(); }
    const Genome *genome(const std::string &name) const;    /* nullptr if absent */

  private:
    std::string _dir;
    Manifest _manifest;
    std::vector<std::unique_ptr<Genome>> _genomes;           /* manifest order */
    std::map<std::string, const Genome *> _byName;
};

/* The lowest common ancestor of two genomes. */
const Genome *lowestCommonAncestor(const Genome *a, const Genome *b);
}

#endif
