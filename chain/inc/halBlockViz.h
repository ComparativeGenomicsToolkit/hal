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
struct block
{
   struct block *next;
   char *qChrom;
   int tStart;
   int qStart;
   int size;
   char strand;
   char *sequence;
};

/** Open a HAL alignment file read-only.  Open the genome specied 
 * by the file for reading.  
 * @param halFileName path to location of HAL file on disk 
 * @param qSpeciesName name of the query species (genome)
 * @param calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return new handle or -1 of open failed.
 *
 * WARNING: NOT THREAD SAFE!!! CLIENT MUST USE LOCKS TO ENSURE THAT
 * THAT THIS METHOD IS NOT CALLED CONCURRENTLY BY DIFFERENT THREADS
*/
int halOpen(char *halFileName, char* qSpeciesName, int threadID);

/** Close a HAL alignment, given its handle
 * @param halHandle previously obtained from halOpen 
 * @param calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return 0: success -1: failure
 *
 * WARNING: NOT THREAD SAFE!!! CLIENT MUST USE LOCKS TO ENSURE THAT
 * THAT THIS METHOD IS NOT CALLED CONCURRENTLY BY DIFFERENT THREADS
 */
  int halClose(int halHandle, int threadID);

/** Free linked list of blocks */
void halFreeBlocks(struct block* block);

/** Create linked list of block structures.  Blocks returned will be all
 * aligned blocks in the query sequence that align to the given range
 * in the reference sequence.  The list will be ordered along the reference.
 * The two immediately adjacent blocks to each aligned query block (adjacent
 * along the query genome) will also be returned if they exist. 
 *
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param tChrom name of the chromosome in reference.
 * @param tStart start position in reference  
 * @param tEnd last + 1 position in reference (if 0, then the size of the 
 * chromosome is used). 
 * @param getSequenceString copy DNA sequence (of query species) into 
 * output blocks if not 0. 
 * @param calling thread ID. No two thread ID's will get a reference
 * to the same alignment instance.  If not using threads just use 0.
 * @return  block structure -- must be freed by halFreeBlocks()
 */
struct block *halGetBlocksInTargetRange(int halHandle,
                                        char* tChrom,
                                        int tStart, int tEnd,
                                        int getSequenceString,
                                        int threadID);
  

#ifdef __cplusplus
}
#endif

#endif
