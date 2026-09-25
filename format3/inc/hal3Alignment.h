/*
 * HAL format 3: reading an alignment.
 *
 * An Alignment is the manifest plus the genome files it names, opened as they are
 * first used, from a pack file or a directory (hal3Store.h), local or remote.  A Genome
 * knows its place in the tree, its sequence table (axis a of its .hgn) and its .hgn,
 * which holds the branch to its parent.  Depends on liftmap only -- no HAL API, no HDF5.
 *
 * Thread safety: an Alignment may be shared by threads; a genome's .hgn is opened
 * once, under a lock, and liftmap files are safe to share.
 */
#ifndef _HAL3ALIGNMENT_H
#define _HAL3ALIGNMENT_H

#include "hal3Manifest.h"
#include "hal3Store.h"
#include "twobit64.h"
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
    mutable std::string _openError;            /* why _hgn could not be opened */
    mutable std::once_flag _dnaOpened;
    mutable lmap_io *_dna = nullptr;           /* the .2bit, for shared openDna handles */
    mutable std::string _dnaError;
    void open() const;
};

class Alignment {
  public:
    /* A format 3 alignment: pack file, directory or manifest, local or a URL.  Throws
     * std::runtime_error. */
    explicit Alignment(const std::string &path);
    ~Alignment();
    Alignment(const Alignment &) = delete;
    Alignment &operator=(const Alignment &) = delete;

    const Manifest &manifest() const { return _manifest; }
    const std::string &path() const { return _path; }
    const Store &store() const { return *_store; }
    /* g's DNA; the caller tb_close()s it, before the Alignment goes.  A shared handle
     * reads through one file handle per genome, opened on first use and kept until the
     * Alignment goes, so reopening is cheap and holds no descriptor of its own (for
     * callers that reopen, like hal maf per chunk); otherwise the handle has its own.
     * Throws if the genome has none. */
    tb_file *openDna(const Genome *g, bool shared = false) const;
    const Genome *root() const { return _genomes.empty() ? nullptr : _genomes[0].get(); }
    const Genome *genome(const std::string &name) const;    /* nullptr if absent */

  private:
    std::string _path;
    std::unique_ptr<Store> _store;
    Manifest _manifest;
    std::vector<std::unique_ptr<Genome>> _genomes;           /* manifest order */
    std::map<std::string, const Genome *> _byName;
};

/* The lowest common ancestor of two genomes. */
const Genome *lowestCommonAncestor(const Genome *a, const Genome *b);
}

#endif
