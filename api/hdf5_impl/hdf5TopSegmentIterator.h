/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5TOPSEGMENTITERATOR_H
#define _HDF5TOPSEGMENTITERATOR_H

#include <H5Cpp.h>
#include "halTopSegmentIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5TopSegment.h"

namespace hal {

class HDF5BottomSegmentIterator;

class HDF5TopSegmentIterator : public TopSegmentIterator
{
   friend class HDF5BottomSegmentIterator;
public:
   
   HDF5TopSegmentIterator(HDF5Genome* genome, hal_index_t index,
                          hal_offset_t startOffset = 0, 
                          hal_offset_t endOffset = 0,
                          hal_bool_t inverted = false);
   ~HDF5TopSegmentIterator();
   
   // ITERATOR METHODS
   void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   void toReverse() const;
   void toSite(hal_index_t position, bool slice = true) const;
   bool hasNextParalogy() const;
   void toNextParalogy() const;
   hal_offset_t getStartOffset() const;
   hal_offset_t getEndOffset() const;
   void slice(hal_offset_t startOffset, hal_offset_t endOffset) const;
   hal_index_t getStartPosition() const;
   hal_size_t getLength() const;
   hal_bool_t getReversed() const;
   void getString(std::string& outString) const;
   bool leftOf(hal_index_t genomePos) const;
   bool rightOf(hal_index_t genomePos) const;
   bool overlaps(hal_index_t genomePos) const;

   //TOP ITERATOR METHODS
   TopSegmentIteratorPtr copy();
   TopSegmentIteratorConstPtr copy() const;
   void copy(TopSegmentIteratorConstPtr ts) const;
   void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const;
   void toParseUp(BottomSegmentIteratorConstPtr bs) const;
   TopSegment* getTopSegment();
   const TopSegment* getTopSegment() const;
   bool equals(TopSegmentIteratorConstPtr other) const;
   bool hasParent() const;
   bool hasParseDown() const;
   bool inRange() const;
protected:
   HDF5TopSegment _topSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable hal_bool_t _reversed;
};

inline bool HDF5TopSegmentIterator::inRange() const
{
  return _topSegment._index >= 0 && _topSegment._index < 
     (hal_index_t)_topSegment._genome->_topArray.getSize();
}

}
#endif
