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

   void addSegment(LodSegment* segment);
   void clear();

   /** Get the total length of all (existing) adjacencies in all segmetns
    * in the block */
   hal_size_t getTotalAdjLength() const;

   /** Extend segments as much as possible (or maxFrac) in both directions. */
   void extend(double maxFrac = 1.0);
      
protected:
  
   /** Get the minimum adjacencly lengths.  
    * NOTE THAT EVERY SEGMENT MUST HAVE AN EXISTING ADJACENCY */
   hal_size_t getMinHeadAdjLen() const;
   hal_size_t getMinTailAdjLen() const;

   SegmentList _segments;

private:
   LodBlock(const LodBlock&);
   const LodBlock& operator=(const LodBlock&) const;
};

inline hal_size_t LodBlock::getNumSegments() const
{
  return _segments.size();
}

inline hal_size_t LodBlock::getLength() const
{
  return _segments.empty() ? 0 : _segments[0]->getLength();
}


}

#endif
