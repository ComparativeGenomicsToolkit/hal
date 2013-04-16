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

/** Blockc struct. 
 * NOTE: ALL COORDINATES ARE FORWARD-STRAND RELATIVE 
 */
struct hal_block_t
{
   struct hal_block_t* next;
   char *qChrom;
   int tStart;
   int qStart;
   int size;
   char strand;
   char *sequence;
};

/** Some information about a genome */
struct hal_species_t
{
   struct hal_species_t* next;
   char* name;
   unsigned length;
   unsigned numChroms;
   char* parentName;
   double parentBranchLength;
};

/** Some information about a sequence */
struct hal_chromosome_t
{
   struct hal_chromosome_t* next;
   char* name;
   unsigned length;
};

/** Open a text file created by halLodInterpolate.py for viewing. 
 * This text file contains a list of paths to progressively coarser
 * levels of detail of a source HAL file.  For example, it may look like
 * 0 ecoli.hal
 * 1000 ecoli_10.hal
 * 10000 ecoli_100.hal
 * This file is saying that for query lengths between 0 and 1000, use 
 * the first (original) file.  For lengths between 1000 and 10000 use 
 * the second file. For lengths greater than 10000 use the third file.
 *
 * halGetBlocksInTargetRange will automatically use the above-described
 * logic.  Calling halOpen (below) is the equivalent of having just one
 * entry (0)
 *
 * @param lodFilePath path to location of HAL LOD file on disk 
 * @return new handle or -1 of open failed.
*/
int halOpenLOD(char *lodFilePath);

/** Open a HAL alignment file read-only.  
 * @param halFilePath path to location of HAL file on disk 
 * @return new handle or -1 of open failed.
*/
int halOpen(char *halFilePath);

/** Close a HAL alignment, given its handle
 * @param halHandle previously obtained from halOpen 
 * @return 0: success -1: failure
 */
int halClose(int halHandle);

/** Free linked list of blocks */
void halFreeBlocks(struct hal_block_t* block);

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
 * @param getSequenceString copy DNA sequence (of query species) into 
 * output blocks if not 0. 
 * @param doDupes create blocks for duplications if not 0.  When this 
 * option is enabled, the same region can appear in more than one block.
 * @return  block structure -- must be freed by halFreeBlocks()
 */
struct hal_block_t *halGetBlocksInTargetRange(int halHandle, 
                                              char* qSpecies,
                                              char* tSpecies,
                                              char* tChrom,
                                              int tStart, int tEnd,
                                              int getSequenceString,
                                              int doDupes);
 

/** Create a linked list of the species in the hal file.
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @return  species structure -- must be freed by client */
struct hal_species_t *halGetSpecies(int halHandle);

/** Create a linked list of the chromosomes in the 
 * @param halHandle handle for the HAL alignment obtained from halOpen 
 * @param speciesName The name of the species whose chromomsomes you want 
 * @return  chromosome structure -- must be freed by client */
struct hal_chromosome_t *halGetChroms(int halHandle, 
                                      char* speciesName);

/** Create a string of the DNA characters of the given range of a chromosome
 * @param halHandle handle for the HAL alignment obtained from halOpen 
 * @param speciesName The name of the species 
 * @param chromName The name of the chromosome within the species
 * @param start The first position of the chromosome
 * @param end The last + 1 position of the chromosome 
 * @return dna string -- must be freed by client */
char *halGetDna(int halHandle,
                char* speciesName,
                char* chromName, 
                int start, int end);

#ifdef __cplusplus
}
#endif

#endif
