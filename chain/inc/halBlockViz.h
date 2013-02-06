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

/** Blockc struct. */
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
 * @return new handle or -1 of open failed.
*/
int halOpen(char *halFileName, char* qSpeciesName);

/** Close a HAL alignment, given its handle
 * @param halHandle previously obtained from halOpen 
 * @return 0: success -1: failure
 */
int halClose(int halHandle);

/** Free linked list of blocks */
void halFreeBlocks(struct block* block);

/** Create linked list of block structures.
 * @param halHandle handle for the HAL alignment obtained from halOpen
 * @param tChrom name of the chromosome in reference.
 * @param tStart start position in reference  
 * @param tEnd last + 1 position in reference  
 * @param copy DNA sequence (of query species) into output blocks
 * @retrun block structure -- must be freed by halFreeBlocks()
 */
struct block *halGetBlocksInTargetRange(int halHandle,
                                        char* tChrom,
                                        int tStart, int tEnd,
                                        int getSequenceString);
  

#ifdef __cplusplus
}
#endif

#endif
