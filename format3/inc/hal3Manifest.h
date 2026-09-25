/*
 * HAL format 3: the manifest.
 *
 * A format 3 alignment is a set of immutable per-genome files plus one small mutable
 * manifest that ties them together: the tree, the genome roster, where each genome's
 * files are and their checksums, and alignment-wide metadata.  Unpacked, it is a
 * directory `foo.hal/` holding `manifest` and `genomes/`; packed (later), the same
 * members in one container file.
 *
 * Each genome has:
 *   - a `.hgn` genome file: a liftmap file (profile "hal.branch") whose axis a is the
 *     genome's sequences and axis b its parent's, with the runs of the branch above it.
 *     The root's has sequences and no runs.
 *   - a `.2bit` with its DNA (optional: DNA may be absent).
 *
 * The manifest is text, one record per line, tab-separated, with \t \n \\ escaped:
 *
 *   #hal-manifest 3
 *   tree    <newick>
 *   meta    <key>    <value>
 *   genome  <name> <parent|-> <branch length|-> <length> <sequences>
 *           <hgn path> <hgn bytes> <hgn crc32> <2bit path|-> <2bit bytes> <2bit crc32>
 *
 * Paths are relative to the manifest's directory.  Lines starting with '#' after the
 * first are comments.  This code depends on nothing in the rest of HAL.
 */
#ifndef _HAL3MANIFEST_H
#define _HAL3MANIFEST_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace hal3 {

extern const char *const MANIFEST_MAGIC;     /* "#hal-manifest 3" */
extern const char *const MANIFEST_NAME;      /* "manifest", inside the alignment directory */
extern const char *const HGN_PROFILE;        /* "hal.branch" */

struct GenomeEntry {
    std::string name;
    std::string parent;          /* "" for the root */
    double branchLength = -1;    /* to the parent; -1 for the root */
    uint64_t length = 0;         /* bases */
    uint64_t numSequences = 0;
    std::string hgnPath;
    uint64_t hgnBytes = 0;
    uint32_t hgnCrc = 0;
    std::string dnaPath;         /* "" when the genome has no DNA */
    uint64_t dnaBytes = 0;
    uint32_t dnaCrc = 0;
};

struct Manifest {
    std::string newick;
    std::map<std::string, std::string> meta;
    std::vector<GenomeEntry> genomes;    /* the root first, then preorder */

    const GenomeEntry *find(const std::string &name) const;

    /* Write to `path` atomically (via path.partial).  Throws std::runtime_error. */
    void write(const std::string &path) const;
    /* Throws std::runtime_error with the line number on malformed input. */
    static Manifest read(const std::string &path);
    /* The same from the manifest's text; what names it in errors. */
    static Manifest parse(const std::string &text, const std::string &what);
    /* The manifest's text, as write() writes it. */
    std::string text() const;
    /* Does the start of a file look like a manifest? */
    static bool isManifest(const std::string &initialBytes);
};

/* A file-system-safe stem for a genome's files: the name, with anything outside
 * [A-Za-z0-9._-] (and a leading '.') percent-encoded. */
std::string genomeFileStem(const std::string &genomeName);

/* CRC-32 (zlib) and size of a whole file.  Throws std::runtime_error. */
uint32_t fileCrc32(const std::string &path, uint64_t *bytes);

/* A parsed Newick tree: each node's name, parent (-1 for the root) and branch
 * length (-1 when absent).  Throws std::runtime_error on malformed input. */
struct NewickNode {
    std::string name;
    int parent;
    double branchLength;
};
std::vector<NewickNode> parseNewick(const std::string &newick);
}

#endif
