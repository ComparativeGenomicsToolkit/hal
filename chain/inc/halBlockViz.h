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

/** Initialize the static memory where C++ pointers to HAL alignments
 * are stored.  We allow 1 pointer per threadID, so maxThreadID alignemnts
 * can be kept
 * @return 0 succes -1 failure */
int halInit(int maxThreadID);

/** Wipe out all the alignment instances kept in the static contained
 * that was created by halInit
 * @return 0 succes -1 failure */
int halExit();

/** Open a HAL alignment file read-only.  
 * @param halFileName path to location of HAL file on disk 
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return new handle or -1 of open failed.
*/
int halOpen(char *halFileName, int threadID);

/** Close a HAL alignment, given its handle
 * @param halHandle previously obtained from halOpen 
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return 0: success -1: failure
 */
int halClose(int halHandle, int threadID);

/** Free linked list of blocks */
void halFreeBlocks(struct hal_block_t* block);

/** Create linked list of block structures.  Blocks returned will be all
 * aligned blocks in the query sequence that align to the given range
 * in the reference sequence.  The list will be ordered along the reference.
 * The two immediately adjacent blocks to each aligned query block (adjacent
 * along the query genome) will also be returned if they exist. 
 *
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
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
                                              int threadID,
                                              char* qSpecies,
                                              char* tSpecies,
                                              char* tChrom,
                                              int tStart, int tEnd,
                                              int getSequenceString,
                                              int doDupes);
 

/** Create a linked list of the species in the hal file.
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return  species structure -- must be freed by client */
struct hal_species_t *halGetSpecies(int halHandle, int threadID);

/** Create a linked list of the chromosomes in the 
 * @param halHandle handle for the HAL alignment obtained from halOpen 
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @param speciesName The name of the species whose chromomsomes you want 
 * @return  chromosome structure -- must be freed by client */
struct hal_chromosome_t *halGetChroms(int halHandle, int threadID,
                                      char* speciesName);

/** Create a string of the DNA characters of the given range of a chromosome
 * @param halHandle handle for the HAL alignment obtained from halOpen 
 * @param threadID calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @param speciesName The name of the species 
 * @param chromName The name of the chromosome within the species
 * @param start The first position of the chromosome
 * @param end The last + 1 position of the chromosome 
 * @return dna string -- must be freed by client */
char *halGetDna(int halHandle, int threadID,
                char* speciesName, char* chromName, 
                int start, int end);

#ifdef __cplusplus
}
#endif

#endif
