/*
 * HAL format 3: a C interface, for programs in C (such as Cactus's cactus_consolidated).
 *
 * Functions return 0 on success and -1 on failure; on failure, if err is given, *err is
 * set to a message the caller frees with free().  Link libHal3.a, twobit64 and liftmap,
 * and the C++ runtime (-lstdc++), zlib and pthreads.
 */
#ifndef _HAL3_C_H
#define _HAL3_C_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- a Cactus subproblem to an alignment fragment (see hal3/hal3Cactus.h) ----------
 *
 * Open, declare every sequence (in .c2h order, with its length), then give each, in the
 * same order, its .c2h segment lines and its DNA; finish.  Genomes the tree does not
 * name (outgroups) are skipped, whatever is given for them. */
typedef struct hal3_c2h hal3_c2h;

/* fragmentDir: a new directory.  newick: the subproblem's ancestor and its ingroup
 * children, e.g. "(a,b)anc;" (branch lengths optional: hal merge --tree can set them).
 * presorted: runs arrive in genome order, as Cactus writes them (streamed, not
 * sorted). */
hal3_c2h *hal3_c2h_open(const char *fragmentDir, const char *newick, int presorted, char **err);
int hal3_c2h_declare(hal3_c2h *h, const char *genome, const char *name, int ancestor, int64_t length, char **err);
int hal3_c2h_sequence(hal3_c2h *h, const char *genome, const char *name, const char *segments, size_t segmentsLength,
                      const char *dna, int64_t dnaLength, char **err);
/* Write the fragment's genome files and manifest.  Frees h, whether or not it succeeds. */
int hal3_c2h_finish(hal3_c2h *h, char **err);
/* Abandon: frees h, leaving no genome file. */
void hal3_c2h_free(hal3_c2h *h);

/* ---- packing ------------------------------------------------------------------------ */

/* The alignment at in (a directory or a pack) as a single-file pack at out. */
int hal3_pack(const char *in, const char *out, char **err);

/* ---- universal columns as trees of bases (see hal3Columns.h) ------------------------
 *
 * The columns of hal maf --universal, duplicates kept, with every base of every genome in
 * scope in one -- ancestors' too -- and each base with its parent: the base of the parent
 * genome it aligns to, which is always in the same column.  A column is so a tree of
 * bases, rooted at its one base with no parent.  The axis is hal maf --axisRow's, so hal
 * wig projects scores given on it.  Open once; then get chunk k (0 <= k < chunks) by
 * worker w (0 <= w < opts.workers): different workers may run at once, from different
 * threads, and chunks may be taken in any order. */
typedef struct hal3_columns hal3_columns;

typedef struct {
    const char *refSequence;      /* only this backbone sequence; NULL: the whole alignment */
    int64_t start, length;        /* with refSequence: length bases from start (length < 0: to its end) */
    const char *rootGenome;       /* the backbone: the subtree of this genome; NULL: the root */
    const char *excludeGenomes;   /* leaves left out, comma-separated; NULL: none */
    int leafInsertions;           /* a leaf's bases with no parent are columns of their own */
    int ancestorBases;            /* 0: ancestors' bases are N (their DNA is not read) */
    int64_t chunk;                /* columns a chunk */
    int workers;                  /* threads that will call hal3_columns_chunk at once */
} hal3_columns_options;

/* the defaults: the whole alignment from the root, ancestor bases, 2^20 columns a chunk,
 * one worker */
void hal3_columns_defaults(hal3_columns_options *opts);

/* A chunk's columns.  The arrays hold one value per entry; column u's entries are
 * [colStart[u], colStart[u+1]), in genome preorder (a parent before its children), the
 * first the column's top.  Valid until worker w's next hal3_columns_chunk, or close. */
typedef struct {
    const char *axisName;      /* the axis sequence the chunk lies on */
    int64_t axisStart;         /* the axis position (0-based) of column 0 */
    int64_t cols;
    const int64_t *colStart;   /* cols + 1 */
    const uint32_t *genome;    /* the base's genome (hal3_columns_genome_name etc.) */
    const int32_t *parent;     /* its parent's entry, counted from colStart[u]; -1: the top */
    const uint32_t *piece;     /* the stretch of its genome it lies in: two entries of a genome under
                                  one parent are copies, and the same pieces over several columns
                                  the same duplication */
    const char *base;          /* in the axis's orientation, soft-masking kept */
    int64_t entries;
} hal3_column_chunk;

hal3_columns *hal3_columns_open(const char *alignment, const hal3_columns_options *opts, char **err);
size_t hal3_columns_chunks(const hal3_columns *c);
int64_t hal3_columns_columns(const hal3_columns *c);   /* in all chunks */
/* the genomes in scope, in preorder from the backbone (index 0) */
uint32_t hal3_columns_genomes(const hal3_columns *c);
const char *hal3_columns_genome_name(const hal3_columns *c, uint32_t i);
int32_t hal3_columns_genome_parent(const hal3_columns *c, uint32_t i);   /* -1 for the backbone */
int hal3_columns_genome_is_leaf(const hal3_columns *c, uint32_t i);
double hal3_columns_branch_length(const hal3_columns *c, uint32_t i);    /* to its parent, from the tree; -1 for the root */
int hal3_columns_chunk(hal3_columns *c, size_t k, int w, hal3_column_chunk *out, char **err);
void hal3_columns_close(hal3_columns *c);

#ifdef __cplusplus
}
#endif

#endif
