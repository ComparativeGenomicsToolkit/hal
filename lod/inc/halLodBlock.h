/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODBLOCK_H
#define _HALLODBLOCK_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <map>
#include "hal.h"
#include "halLodSegment.h"

namespace hal {

class LodBlock;

std::ostream& operator<<(std::ostream& os, const LodBlock& block);

struct LodBlockPBigger
{
   bool operator()(const LodBlock* b1, const LodBlock* b2) const;
};

/* A block is a list of homolgous segments.  All these segments must
 * be the same length.  The block's destructor will free all segments
 * it contains. 
 */
class LodBlock
{
   friend std::ostream& operator<<(std::ostream& os, const LodBlock& block);

public:

   typedef std::vector<LodSegment*> SegmentList;
   typedef SegmentList::iterator SegmentIterator;
   typedef SegmentList::const_iterator SegmentConstIterator;

   LodBlock();
   ~LodBlock();

   hal_size_t getNumSegments() const;
   hal_size_t getLength() const;
   const LodSegment* getSegment(hal_size_t index) const;

   void addSegment(LodSegment* segment);
   void clear();

   /** Get the total length of all (existing) adjacencies in all segmetns
    * in the block */
   hal_size_t getTotalAdjLength() const;

   /** Extend segments as much as possible (or maxFrac) in both directions. */
   void extend(double maxFrac = 1.0);

   /** Test if all segments in block have head to tail adjacencies
    * to segments in the same non-telomere block, and that all these 
    * adjacencies have length 0.  If test passes, return the candidate
    * adjacent block.  NULL otherwise */
   LodBlock* getHeadMergePartner();

   /** Merge head of this block to tail of adjBlock (which was found with
    * getHeadMergePartner.  Merged segments will disappear and need
    * to be accounted for elsewhere */
   void mergeHead(LodBlock* adjBlock);
   
   /** Insert new blocks as neighbours until all adjacencies have length
    * 0.  (if there are no self edges, at most 1 head block and 1 tail
    * block are created.  If there are self edges, it can take multiple
    * blocks to reduce all the edge  lengths */
   void insertNeighbours(std::vector<LodBlock*>& outList);
      
protected:

   /** Create a new block and insert it as a neighbour.  All adjacencies
    * of this block become 0. */
   LodBlock* insertNewTailNeighbour();
   LodBlock* insertNewHeadNeighbour();
  
   /** Get the maximum length to extend the block.  This is equivalent
    * to the minimum adjacency length, except that adjacencies between
    * two segments within this blocks are cut in half (since extending
    * eats them in both directions. 
    * NOTE THAT EVERY SEGMENT MUST HAVE AN EXISTING ADJACENCY */
   hal_size_t getMaxHeadExtensionLen() const;
   hal_size_t getMaxTailExtensionLen() const;

   /** Get the maximum lengths for insertion.  This is the equivalent
    * to the minimum non-zero adjacency.  If no nonzero lengths exist
    * then 0 is returned */
   hal_size_t getMaxHeadInsertionLen() const;
   hal_size_t getMaxTailInsertionLen() const;

   SegmentList _segments;

private:
   LodBlock(const LodBlock&);
   const LodBlock& operator=(const LodBlock&) const;
};

inline bool LodBlockPBigger::operator()(const LodBlock* b1, 
                                        const LodBlock* b2) const
{
  bool bigger = b1->getNumSegments() > b2->getNumSegments();
  return bigger;
}

inline hal_size_t LodBlock::getNumSegments() const
{
  return _segments.size();
}

inline hal_size_t LodBlock::getLength() const
{
  return _segments.empty() ? 0 : _segments[0]->getLength();
}

inline const LodSegment* LodBlock::getSegment(hal_size_t index) const
{
  assert(index < getNumSegments());
  return _segments.at(index);
}

}

#endif
