/*
 * HAL format 3: the manifest.
 *
 * A format 3 alignment is a set of immutable per-genome files plus one small mutable
 * manifest that ties them together: the tree, the genome roster, where each genome's
 * files are and their checksums, and alignment-wide metadata.  Unpacked, it is a
 * directory `foo.hal/` holding `manifest` and `genomes/`; packed, the same members in
 * one file (hal3Store.h).
 *
 * Each genome has:
 *   - a `.hgn` genome file: a liftmap file (profile "hal.branch") whose axis a is the
 *     genome's sequences and axis b its parent's, with the runs of the branch above it.
 *     The root's has sequences and no runs.
 *   - a `.2bit` with its DNA (optional: DNA may be absent).
 *   - optionally, other files about its own sequences ("members"), each with a role:
 *     "lik", a `.hlk` of per-base likelihoods (hal3Likelihood.h).  Tools that do not
 *     know a role carry its file unchanged.
 *
 * The manifest is text, one record per line, tab-separated, with \t \n \\ escaped:
 *
 *   #hal-manifest 3
 *   tree    <newick>
 *   meta    <key>    <value>
 *   genome  <name> <parent|-> <branch length|-> <length> <sequences>
 *           <hgn path> <hgn bytes> <hgn crc32> <2bit path|-> <2bit bytes> <2bit crc32>
 *   member  <genome> <role> <path> <bytes> <crc32>     (after its genome's line)
 *
 * Paths are member names, `genomes/<name>` (no other '/', no "." or ".."), each named
 * once: relative to the manifest's directory, or in a pack the names in its index.  Lines
 * starting with '#' after the first are comments.  This code depends on nothing in the
 * rest of HAL.
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

/* A genome's file beyond its .hgn and .2bit. */
struct MemberRef {
    std::string role;            /* e.g. "lik" */
    std::string path;
    uint64_t bytes = 0;
    uint32_t crc = 0;
};

/* The file name extension a member of this role gets: ".hlk" for "lik", else "." + role. */
std::string memberExtension(const std::string &role);
/* A role is lower-case letters, digits and '_', starting with a letter, and never
 * "hgn", "2bit", "hlk" or "dna" -- so a member's file never has a genome's own file's
 * name, nor another role's. */
bool validMemberRole(const std::string &role);

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
    std::vector<MemberRef> members;

    const MemberRef *member(const std::string &role) const;   /* nullptr if none */
};

/* A file the manifest names, with its size and CRC-32. */
struct NamedFile {
    std::string path;
    uint64_t bytes;
    uint32_t crc;
};

struct Manifest {
    std::string newick;
    std::map<std::string, std::string> meta;
    std::vector<GenomeEntry> genomes;    /* the root first, then preorder */

    const GenomeEntry *find(const std::string &name) const;
    /* Every file the manifest names -- each genome's .hgn, .2bit and members -- in
     * genome order (not the manifest itself). */
    std::vector<NamedFile> files() const;

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
