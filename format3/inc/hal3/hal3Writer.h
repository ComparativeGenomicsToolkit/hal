/*
 * HAL format 3 API: writing alignments.
 *
 * An alignment is written genome by genome, independently: a GenomeWriter makes one
 * genome's .hgn (its sequences, its parent's, and the runs of the branch between) and
 * its .2bit, checking as it goes.  A Fragment gathers such genomes into a small
 * alignment -- for Cactus, one subproblem: an ancestor as its root and its children.
 * merge() stitches fragments whose trees meet in shared genomes (a subproblem's root is
 * a leaf of its parent's) into one alignment, writing each file once.  validate() is
 * hal validate's check, as a call.
 *
 * Nothing here needs another genome's files: the runs of a branch name the parent's
 * sequences by index, so a writer needs the parent's sequence table (names, lengths),
 * not the parent.  Files are written as <name>.partial and renamed at the end, so a
 * killed job never leaves a half-written file under a final name.
 */
#ifndef _HAL3_API_WRITER_H
#define _HAL3_API_WRITER_H

#include "hal3Error.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace hal3 {
namespace api {

struct SequenceDef {
    std::string name;
    int64_t length;
};

struct GenomeWriterOptions {
    uint32_t chunkRuns = 8192;           /* runs per chunk at most */
    uint64_t chunkParentSpan = 1000000;  /* parent bases a chunk spans at most */
    /* runs promised in child order -- each child sequence's together, starts
     * increasing -- so they stream to the file with nothing held or sorted; a run that
     * breaks the promise fails */
    bool presorted = false;
    size_t memBytes = 0;                 /* runs held before spilling to disk (0: 1 GiB) */
    std::string tmpDir;                  /* where they spill ("": beside the output) */
    std::string source = "hal3";         /* the file's "source" metadata: what wrote it */
    bool checksums = true;               /* finish() computes each file's CRC-32 */
};

/* What a finished genome is, as a manifest records it. */
struct GenomeFiles {
    std::string name, parent;            /* parent "" for a root */
    std::string hgnFile, dnaFile;        /* the files written; dnaFile "" without DNA */
    int64_t length = 0;
    uint64_t numSequences = 0, runs = 0;
    uint64_t hgnBytes = 0, dnaBytes = 0;
    uint32_t hgnCrc = 0, dnaCrc = 0;     /* 0 unless checksums */
};

class GenomeWriter {
  public:
    /* Writes dir/<stem>.hgn (and dir/<stem>.2bit if DNA is added), stem being the
     * genome's file stem (the name, made file-system safe).  parent "" for a root, which
     * has no parent sequences and no runs.  Throws Error(Argument) on duplicate sequence
     * names, Error(Io) if the files cannot be written. */
    GenomeWriter(const std::string &dir, const std::string &name, const std::string &parent,
                 const std::vector<SequenceDef> &sequences, const std::vector<SequenceDef> &parentSequences,
                 const GenomeWriterOptions &opts = GenomeWriterOptions());
    /* Unfinished: everything is discarded, nothing left behind. */
    ~GenomeWriter();
    GenomeWriter(const GenomeWriter &) = delete;
    GenomeWriter &operator=(const GenomeWriter &) = delete;

    /* the genome's own metadata (HAL's per-genome metadata), any time before finish */
    void setMeta(const std::string &key, const std::string &value);
    /* [childStart, +length) of child sequence childSeq aligned to [parentStart, +length)
     * of parent sequence parentSeq; reversed pairs childStart+i with
     * parentStart+length-1-i.  Runs may come in any order (unless presorted), and may
     * split or repeat a colinear stretch: they are merged.  A child base aligned to two
     * parent bases fails finish().  Throws Error(Argument) for a run outside its
     * sequences. */
    void addRun(uint32_t childSeq, int64_t childStart, uint32_t parentSeq, int64_t parentStart, int64_t length, bool reversed);
    /* DNA, every sequence in order, each exactly its declared length: whole, or
     * piecewise.  A to Z in either case (lower case: soft-masked); letters other than
     * ACGT are stored as N.  None at all: no .2bit. */
    void addDna(const char *bases, int64_t length);
    void beginDna();
    void appendDna(const char *bases, size_t n);
    void endDna();
    /* Close both files.  Throws Error(Format) if the runs give a child base two parents
     * or not every sequence got its DNA; nothing is left behind then. */
    GenomeFiles finish();

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

/* A small alignment written genome by genome, in a new directory. */
class Fragment {
  public:
    /* Creates dir and dir/genomes; dir must not already hold an alignment. */
    explicit Fragment(const std::string &dir);
    ~Fragment();
    Fragment(const Fragment &) = delete;
    Fragment &operator=(const Fragment &) = delete;

    /* where GenomeWriters for this fragment write */
    std::string genomesDir() const;
    /* a finished genome, with the length of the branch above it (ignored for the root);
     * the tree follows the genomes' parents, children in the order added */
    void add(const GenomeFiles &genome, double branchLength = -1);
    void setMeta(const std::string &key, const std::string &value);
    /* Write the manifest (one root; every other genome's parent in the fragment), then
     * check the fragment as validate() does, without checksums.  Throws Error(Format)
     * with the problems. */
    void finish();

  private:
    struct Impl;
    std::unique_ptr<Impl> _impl;
};

struct MergeOptions {
    bool pack = true;                    /* out a single-file pack; else a directory */
    bool link = true;                    /* directory out: hard-link files from directory
                                          * fragments where possible, else copy */
    bool checkDna = true;                /* a genome's DNA in two fragments must agree */
    bool checksums = false;              /* validate the result with checksums too */
};

/* Stitch fragments (packs or directories) into one alignment at out, which must not
 * exist.  Their trees must meet only in shared genomes, into one tree with one root.  A
 * genome in several fragments takes its .hgn from the one where it has a parent, its
 * .2bit from the one where it is the root (else any): each checked against the
 * fragment's manifest as it is copied.  The tree splices the fragments' at shared
 * genomes, children in fragment order.  Metadata: the fragments', which must not
 * conflict, then meta over them.  The result is validated; on any problem nothing is
 * left at out and Error(Format) says why. */
void merge(const std::vector<std::string> &fragments, const std::string &out,
           const std::map<std::string, std::string> &meta = std::map<std::string, std::string>(),
           const MergeOptions &opts = MergeOptions());

struct ValidateOptions {
    bool checksums = true;               /* every file's CRC-32 against the manifest */
    bool deep = true;                    /* liftmap's full check of every .hgn, .2bit and
                                          * likelihood chunk; else the tables only */
};

/* Every problem with the alignment at path (a pack, directory or manifest), empty when
 * it is valid: the manifest and tree; each genome's files, their names, sizes and
 * checksums; axis a a function (one parent per base); each genome's sequences, its
 * .2bit's and its children's view of them the same; no stray members. */
std::vector<std::string> validate(const std::string &path, const ValidateOptions &opts = ValidateOptions());
}
}

#endif
