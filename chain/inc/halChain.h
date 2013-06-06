/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCHAIN_H
#define _HALCHAIN_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"

namespace hal {

// Chain structure as described in 
//  http://genome.ucsc.edu/goldenPath/help/chain.html
// defines pairwise alignment as list colinear blocks. 
struct ChainBlock
{
   hal_size_t _size;
   hal_size_t _tGap;
   hal_size_t _qGap;
};

struct Chain
{
   std::string _tName;
   hal_size_t _tSize;
   hal_size_t _tStart;
   hal_size_t _tEnd;
   char _tStrand;

   std::string _qName;
   hal_size_t _qSize;
   hal_size_t _qStart;
   hal_size_t _qEnd;
   char _qStrand;
   hal_size_t _id;
   std::vector<ChainBlock> _blocks;
};

std::ostream& operator<<(std::ostream&, const ChainBlock& b);
std::ostream& operator<<(std::ostream&, const Chain& c);

/** Convert a gapped iterator to a chain (with resepct to its parent in the
 * hal graph).   
 * @param gt Gapped Top iterator to convert
 * @param outChain chain structure to write (or overwrite)
 * @param startOffset bases to skip at beginning (gapped iterators don't 
 * support this natively.
 * @param endOffset bases to skip at end (gapped iterators don't support this
 * natively.
 */
void gtIteratorToChain(GappedTopSegmentIteratorConstPtr gt, Chain& outChain,
                       hal_offset_t startOffset = 0, 
                       hal_offset_t endOffset = 0);

}



#endif
