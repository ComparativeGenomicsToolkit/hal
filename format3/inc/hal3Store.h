/*
 * HAL format 3: where an alignment's members live.
 *
 * An alignment is a set of named members -- "manifest", and "genomes/<stem>.hgn" and
 * "genomes/<stem>.2bit" per genome -- stored either as
 *
 *   a directory:  the members are files under it (the development and working form), or
 *   a pack:       one file holding them all (the form to use and share).
 *
 * Both open the same way, locally or over HTTP(S), and every member is read through
 * liftmap's lmap_io, so the genome files never know which it is.
 *
 * The pack format:
 *
 *   [0, 512)   header: magic "\x89HAL\r\n\x1a\n", u32 format 3, u32 pack version 1, then
 *              two commit slots of 32 bytes each, little-endian:
 *                u64 generation, u64 index offset, u64 index length,
 *                u32 CRC-32 of the index, u32 CRC-32 of the slot's first 28 bytes
 *   members    their bytes, anywhere after the header
 *   index      text: "#hal3-pack-index 1", then one line per member,
 *              name <tab> offset <tab> length <tab> CRC-32 (8 hex digits)
 *
 * The current version is the valid slot with the higher generation (whose index checks
 * out).  A change appends -- new members, then a new index -- flushes that to disk, and
 * only then overwrites the older slot, and flushes again: a crash at any point leaves
 * the previous version intact and current.  Members a change replaces or removes stay
 * in the file as dead space until the pack is rewritten (hal pack).  Nothing is ever
 * written into a member a reader may be using.
 */
#ifndef _HAL3STORE_H
#define _HAL3STORE_H

#include "lmap_io.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace hal3 {

extern const char PACK_MAGIC[8];

class Store {
  public:
    /* A directory, a pack file, or a manifest file (its directory); local or a URL.
     * Throws std::runtime_error. */
    static std::unique_ptr<Store> open(const std::string &path);
    virtual ~Store() {}

    virtual bool isPack() const = 0;
    const std::string &path() const { return _path; }

    /* Member names.  For a directory, what is there (for a remote directory, only
     * what the manifest names: a URL cannot be listed). */
    virtual std::vector<std::string> members() const = 0;
    virtual bool has(const std::string &name) const = 0;
    virtual int64_t size(const std::string &name) const = 0;
    /* A handle on one member; the caller closes it, before the store goes. */
    virtual lmap_io *openMember(const std::string &name) const = 0;

    /* The whole of a (small) member. */
    std::string read(const std::string &name) const;
    /* CRC-32 (zlib) of a member's bytes: reads it all. */
    uint32_t crc(const std::string &name) const;

  protected:
    std::string _path;
};

/* A pack's committed state, for listing. */
struct PackEntry {
    uint64_t offset, length;
    uint32_t crc;
};
struct PackInfo {
    uint64_t generation = 0, fileBytes = 0, indexOffset = 0, indexLength = 0;
    std::map<std::string, PackEntry> entries;
};
/* The committed state of a pack store; throws if the store is not a pack. */
const PackInfo &packInfo(const Store &store);

/* Write the alignment at in (a directory or a pack) as a new pack at out: the members its
 * manifest names, each checked against the manifest's CRC-32 as it is copied.  out may be
 * in.  Returns the bytes copied.  Throws std::runtime_error, leaving out as it was. */
uint64_t packAlignment(const std::string &in, const std::string &out);

/* Writing a pack: members are added (appended) and removed, then commit() makes the
 * change current, atomically.  One writer at a time (an exclusive lock on the file). */
class PackWriter {
  public:
    /* A new pack, written as path.partial and renamed into place by the first commit. */
    static std::unique_ptr<PackWriter> create(const std::string &path);
    /* An existing pack, to change. */
    static std::unique_ptr<PackWriter> update(const std::string &path);
    ~PackWriter();

    /* Add or replace a member: from a file, from bytes, or from any lmap_io (such as
     * another store's member).  The data goes to disk now; it is part of the pack only
     * once committed.  Returns the member's CRC-32. */
    uint32_t putFile(const std::string &name, const std::string &path);
    uint32_t putBytes(const std::string &name, const std::string &bytes);
    uint32_t putIo(const std::string &name, lmap_io *io);
    void remove(const std::string &name);
    /* to names the bytes that from named when the writer was opened -- whatever has
     * been done to either name since -- where they are: nothing is copied.  from is
     * left as it is; remove() it if it is to go.  (So a cycle of renames, a swap, can
     * be done one link at a time.) */
    void link(const std::string &from, const std::string &to);
    bool has(const std::string &name) const { return _entries.count(name) > 0; }
    std::vector<std::string> members() const;

    /* Write the index and switch to it. */
    void commit();

  private:
    PackWriter() {}
    std::string _path, _writePath;
    int _fd = -1;
    bool _renamePending = false;
    uint64_t _end = 0, _generation = 0;
    int _nextSlot = 0;
    std::map<std::string, PackEntry> _entries;
    std::map<std::string, PackEntry> _opened;   /* _entries as the writer found them */
    void append(const void *data, size_t n, uint64_t at);
};
}

#endif
