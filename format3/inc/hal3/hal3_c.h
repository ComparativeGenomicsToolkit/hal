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

#ifdef __cplusplus
}
#endif

#endif
