/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBLOCKVIZ_H
#define _HALBLOCKVIZ_H

#ifdef __cplusplus
extern "C" {
#endif

/** This is all prototype code to evaluate how to get blocks streamed 
 * from HAL to the browser. Interface is speficied by Brian */

/* keep integer type definition in one place */
typedef long hal_int_t;

/** range of coordinates in target */
struct hal_target_range_t
{
   struct hal_target_range_t* next;
   hal_int_t tStart;
   hal_int_t size;
};

/** paralgous ranges in the *target* genome, which can't be displayed
 * as snakes so they get filtered out into a separate API */
struct hal_target_dupe_list_t
{
   struct hal_target_dupe_list_t* next;
   hal_int_t id;  
   struct hal_target_range_t* tRange;
   char* qChrom;
};

/** Contains mapped blocks along with target paralgous "blue line" blocks
 * in a separate list.  No overlapping query blocks will be present in
 * mappedBlocks.  If there are three paralgous reference segemnts, then
 * there will be three instances in targetDupeBlocks, but only one
 * "cannonical" instance in mappedBlocks */
struct hal_block_results_t
{
   struct hal_block_t* mappedBlocks;
   struct hal_target_dupe_list_t* targetDupeBlocks;
};
  
/** Blockc struct. 
 * NOTE: ALL COORDINATES ARE FORWARD-STRAND RELATIVE 
 */
struct hal_block_t
{
   struct hal_block_t* next;
   char *qChrom;
   hal_int_t tStart;
   hal_int_t qStart;
   hal_int_t size;
   char strand;
   char *qSequence; // query DNA, if requested
   char *tSequence; // target DNA, if requested
};

/** Some information about a genome */
struct hal_species_t
{
   struct hal_species_t* next;
   char* name;
   hal_int_t length;
   hal_int_t numChroms;
   char* parentName;
   double parentBranchLength;
};

/** Some information about a sequence */
struct hal_chromosome_t
{
   struct hal_chromosome_t* next;
   char* name;
   hal_int_t length;
};

/** Represents genome or alignment metadata as a linked list. */
struct hal_metadata_t
{
  struct hal_metadata_t *next;
  char *key;
  char *value;
};

/** Duplication mode toggler.  
 * HAL_NO_DUPS: No duplications computed
 * HAL_QUERY_DUPS: The same query range can map to multiple places in target
 * HAL_QUERY_AND_TARGET_DUPS: As above, but target can also map to multiple
 * places in query.  Note that these instances are returned in the special
 * target dup list structure */
typedef enum
{
  HAL_NO_DUPS = 0,
  HAL_QUERY_DUPS,
  HAL_QUERY_AND_TARGET_DUPS,
} hal_dup_type_t;

/** Sequence string mode toggler.  
 * HAL_NO_SEQUENCE: Do not return DNA sequence
 * HAL_LOD0_SEQUENCE: Only return sequence if query maps to Lod-0
 * HAL_FORCE_LOD0_SEQUENCE: Return DNA sequence no matter what by forcing
 * query to LOD0 */
typedef enum
{
  HAL_NO_SEQUENCE = 0,
  HAL_LOD0_SEQUENCE,
  HAL_FORCE_LOD0_SEQUENCE,
} hal_seqmode_type_t;

/** Open a text file created by halLodInterpolate.py for viewing. 
 * This text file contains a list of paths to progressively coarser
 * levels of detail of a source HAL file.  For example, it may look like
 * 0 ecoli.hal
 * 1000 lod/ecoli_10.hal
 * 10000 lod/ecoli_100.hal
 * This file is saying that for query lengths between 0 and 1000, use 
 * the first (original) file.  For lengths between 1000 and 10000 use 
 * the second file. For lengths greater than 10000 use the third file.
 *
 * NOTE: If the hal file paths are relative (do not begin with /) as 
 * they are in the above example, then they are assumed to be relative
 * to the directory containing lodFilePath.  If they are absolute, then
 * they will be read as-is. Paths that contain ":/" are assumed to be
 * web addressed of some sort and considered absolute.
 *
 * halGetBlocksInTargetRange will automatically use the above-described
 * logic.  Calling halOpen (below) is the equivalent of having just one
 * entry (0)
 *
 * @param lodFilePath path to location of HAL LOD file on disk 
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return new handle or -1 of open failed.
*/
int halOpenLOD(char *lodFilePath, char **errStr);

/** Open a HAL alignment file read-only.  
 * @param halFilePath path to location of HAL file on disk 
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return new handle or -1 of open failed.
*/
int halOpen(char *halFilePath, char **errStr);

/** Close a HAL alignment, given its handle
 * @param halHandle previously obtained from halOpen
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return 0: success -1: failure
 */
int halClose(int halHandle, char **errStr);

/** Free block results structure */
void halFreeBlockResults(struct hal_block_results_t* results);

/** Free linked list of blocks */
void halFreeBlocks(struct hal_block_t* block);

/** Free linked list of dupe lists*/
void halFreeTargetDupeLists(struct hal_target_dupe_list_t* dupes);

/** Free linked list of species */
void halFreeSpeciesList(struct hal_species_t *species);

/** Get the list of species that can act as coalescence limits when
 * mapping from qSpecies to tSpecies, starting from the MRCA and
 * traveling up.
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param qSpecies the name of the query species.
 * @param tSpecies the name of the reference species.
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return species list -- must be freed by halFreeSpeciesList().
 * NULL on failure.
*/
struct hal_species_t *halGetPossibleCoalescenceLimits(int halHandle,
                                                      const char *qSpecies,
                                                      const char *tSpecies,
                                                      char **errStr);

/** Create linked list of block structures.  Blocks returned will be all
 * aligned blocks in the query sequence that align to the given range
 * in the reference sequence.  The list will be ordered along the reference.
 * The two immediately adjacent blocks to each aligned query block (adjacent
 * along the query genome) will also be returned if they exist. 
 *
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param qSpecies the name of the query species.
 * @param tSpecies the name of the reference species.
 * @param tChrom name of the chromosome in reference.
 * @param tStart start position in reference  
 * @param tEnd last + 1 position in reference (if 0, then the size of the 
 * chromosome is used). 
 * @param tReversed Input region is on the reverse strand (but still 
 * in forward coordinates like BED).  can only be used in liftOverMode
 * otherwise must be set to 0
 * @param seqMode Specify logic used to determine whether or not DNA 
 * sequence is returned (details in hal_seqmode_type_t comments)  
 * @param dupMode Specifies which types of duplications to compute. 
 * (note that when tReversed != 0, target duplications are not supported,
 * so when doing liftover use no dupes or query dupes only) 
 * @param mapBackAdjacencies Species if segments adjacent to query segments
 * in the alignment are mapped back ot the target (producing offscreen 
 * results). Must be set to 0 it tReversed is not 0 (so set to 0 when doing 
 * liftover).
 * @param coalescenceLimitName The name of the genome to use as the
 * limit for walking back and capturing additional paralogs. If NULL,
 * defaults to using the MRCA, which gives no extra paralogs.
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return  block structure -- must be freed by halFreeBlockResults().
 * NULL on failure.
 */
struct hal_block_results_t *halGetBlocksInTargetRange(int halHandle, 
                                                      char* qSpecies,
                                                      char* tSpecies,
                                                      char* tChrom,
                                                      hal_int_t tStart, 
                                                      hal_int_t tEnd,
                                                      hal_int_t tReversed,
                                                      hal_seqmode_type_t seqMode,
                                                      hal_dup_type_t dupMode,
                                                      int mapBackAdjacencies,
                                                      const char *coalescenceLimitName,
                                                      char **errStr);

/** Read alignment into an output file in MAF format.  Interface very 
 * similar to halGetBlocksInTargetRange except multiple query species 
 * can be specified
 *
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param qSpeciesNames the names of the query species (no other information
 * is read from the hal_species_t structure -- just the names).
 * @param tSpecies the name of the reference species.
 * @param tChrom name of the chromosome in reference.
 * @param tStart start position in reference  
 * @param tEnd last + 1 position in reference (if 0, then the size of the 
 * chromosome is used). 
 * @param doDupes create blocks for duplications if not 0.  When this 
 * option is enabled, the same region can appear in more than one block.
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return  number of bytes written. -1 on failure.
 */
hal_int_t halGetMAF(FILE* outFile,
                    int halHandle, 
                    struct hal_species_t* qSpeciesNames,
                    char* tSpecies,
                    char* tChrom,
                    hal_int_t tStart, 
                    hal_int_t tEnd,
                    int doDupes,
                    char **errStr);

/** Create a linked list of the species in the hal file.
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return  species structure -- must be freed by client. NULL on failure. */
struct hal_species_t *halGetSpecies(int halHandle, char **errStr);

/** Create a linked list of the chromosomes in the
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param speciesName The name of the species whose chromomsomes you want
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return chromosome structure -- must be freed by client. NULL on failure */
struct hal_chromosome_t *halGetChroms(int halHandle,
                                      char* speciesName,
                                      char **errStr);

/** Create a string of the DNA characters of the given range of a chromosome
 * @param halHandle handle for the HAL alignment obtained from halOpen 
 * @param speciesName The name of the species 
 * @param chromName The name of the chromosome within the species
 * @param start The first position of the chromosome
 * @param end The last + 1 position of the chromosome 
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return dna string -- must be freed by client. NULL on failure. */
char *halGetDna(int halHandle,
                char* speciesName,
                char* chromName, 
                hal_int_t start, hal_int_t end,
                char **errStr);

/** Get the maximum query size supported by the lod.txt file.  Queries > than
 * this length will return an error. 
 * @param  halHandle handle for the HAL LOD.txt obtained from halOpenLOD 
 * @param errStr pointer to a string that contains an error message on
 * failure. If NULL, throws an exception on failure instead.
 * @return Maximum query length.  Any query > than this value will be invalid.
 *         In the event of an error, -1 will be returned. */
hal_int_t halGetMaxLODQueryLength(int halHandle, char **errStr);

/** Get the metadata for the genome as a linked list instead of a hash.
    Returns NULL if there isn't any metadata for this genome. */
struct hal_metadata_t *halGetGenomeMetadata(int halHandle,
                                            const char *genomeName,
                                            char **errStr);

/** Free a metadata linked list. */
void halFreeMetadataList(struct hal_metadata_t *metadata);

/** Free a linked list of chromosome info. */
void halFreeChromList(struct hal_chromosome_t *chromosome);

#ifdef __cplusplus
}
#endif

#endif
