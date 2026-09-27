/*
 * HAL format 3: changing an alignment's genomes -- splitting, joining, subsetting,
 * replacing, renaming -- from alignments that already exist.  (New mappings come from
 * aligning, outside HAL; this combines and edits what is there.)
 *
 * Operations act on a working model of the tree -- which genomes, their parents,
 * branch lengths, and where each genome's files come from (this alignment or a donor).
 * They apply in the order given and may pass through states that are not valid (a
 * genome removed before the genome that replaces its children arrives); check() then
 * requires the result to be a valid alignment, and write() / writeInPlace() refuse
 * otherwise:
 *
 *   one root, every other genome's parent present;
 *   each genome's .hgn holds its parent's sequence table, name for name and length
 *   for length, as its axis b (the root's none).
 *
 * Writing: every genome's file is compared with the genome's final place in the tree
 * (its name, its parent's name, whether it has a branch).  A file that does not match
 * -- a renamed genome, or its children, or a genome that became the root -- is rewritten;
 * every other file is copied untouched (or, in place, not touched at all).
 *
 * A genome's members (its likelihoods, say: hal3Manifest.h) are about its own bases, so
 * they go wherever its DNA goes: kept, renamed or dropped with it, taken from a donor
 * with its DNA.  They are never rewritten.  Likelihoods that an update makes out of date
 * -- a genome's "subtree" likelihoods when anything below it changes, "posterior" ones
 * when the tree or any branch does -- are kept, and noted in the log.
 */
#ifndef _HAL3UPDATE_H
#define _HAL3UPDATE_H

#include "hal3Manifest.h"
#include "hal3Store.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace hal3 {

class Update {
  public:
    /* The alignment to change (pack or directory).  Throws std::runtime_error. */
    explicit Update(const std::string &path);
    ~Update();

    /* Remove genomes.  A genome with children leaves them without a parent, which the
     * same update must resolve (remove, graft or replace them). */
    void remove(const std::vector<std::string> &genomes);
    /* Remove a genome and everything below it. */
    void prune(const std::string &genome);
    /* Keep only the subtree of genome, which becomes the root. */
    void extract(const std::string &genome);
    /* Keep only these genomes and the ancestors that connect them (their spanning tree),
     * rooted at their common ancestor. */
    void keep(const std::vector<std::string> &genomes);
    /* Join: the donor alignment's root must be a genome here; the donor's other genomes
     * are added below it. */
    void graft(const std::string &donorPath);
    /* Take these genomes' files (and parent, and branch length) from the donor. */
    void replace(const std::string &donorPath, const std::vector<std::string> &genomes);
    void rename(const std::string &from, const std::string &to);
    void setBranchLength(const std::string &genome, double length);
    /* A genome's member of this role (e.g. "lik") from a local file, replacing any it
     * has; or none. */
    void setMember(const std::string &genome, const std::string &role, const std::string &file);
    void dropMember(const std::string &genome, const std::string &role);
    /* A genome's DNA from a local .2bit file (the same sequences). */
    void setDna(const std::string &genome, const std::string &file);

    /* Every problem with the result, empty if it is a valid alignment (and then any
     * likelihoods it makes out of date are noted in the log). */
    std::vector<std::string> check();
    /* A new alignment at out (a pack, or with directory a directory).  Throws, writing
     * nothing, if check() finds problems. */
    void write(const std::string &out, bool directory);
    /* Change the input itself: for a pack, one crash-safe transaction; for a directory,
     * new files, then the manifest, then removals. */
    void writeInPlace();

    /* What the update changes, one line each. */
    const std::vector<std::string> &log() const { return _log; }

  private:
    struct Source {
        const Store *store;             /* the input or a donor; nullptr: a local file */
        std::string member;             /* "" when there is no such file (no DNA) */
    };
    struct Node {
        std::string name, parent;
        double branchLength = -1;
        std::vector<std::string> children;
        Source hgn, dna;
        std::map<std::string, Source> members;   /* by role */
    };
    std::string _path;
    std::unique_ptr<Store> _in;
    std::vector<std::unique_ptr<Store>> _donors;
    Manifest _manifest;             /* the input's, for its metadata and tree text */
    std::map<std::string, Node> _nodes;
    std::string _root;
    bool _treeChanged = false;
    bool _branchesChanged = false;          /* anything a posterior depends on */
    std::set<std::string> _changedBelow;    /* genomes with something changed below them */
    std::vector<std::string> _log;

    static lmap_io *openSource(const Source &s);   /* a store's member, or a local file */
    static bool sourceExists(const Source &s);
    Node &need(const std::string &name);
    void changedBelow(const std::string &genome);   /* genome and its ancestors */
    void staleNotes();
    void detach(const std::string &name);
    void dropSubtree(const std::string &name);
    std::vector<std::string> preorder() const;
    std::string newick() const;
    Manifest finalManifest() const;
};
}

#endif
