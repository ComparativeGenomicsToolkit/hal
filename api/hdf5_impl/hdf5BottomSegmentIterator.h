/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5BOTTOMSEGMENTITERATOR_H
#define _HDF5BOTTOMSEGMENTITERATOR_H

#include <limits>
#include <H5Cpp.h>
#include "halBottomSegmentIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5BottomSegment.h"

namespace hal {

class HDF5TopSegmentIterator;

class HDF5BottomSegmentIterator : public BottomSegmentIterator
{

   friend class HDF5TopSegmentIterator;

public:
   
   HDF5BottomSegmentIterator(HDF5Genome* genome, hal_index_t index,
                             hal_size_t startOffset = 0, 
                             hal_size_t endOffset = 0,
                             hal_bool_t inverted = false);
   ~HDF5BottomSegmentIterator();
   
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

   //BOTTOM ITERATOR METHODS
   BottomSegmentIteratorPtr copy();
   BottomSegmentIteratorConstPtr copy() const;
   void copy(BottomSegmentIteratorConstPtr bs) const;
   void toParent(TopSegmentIteratorConstPtr ts) const; 
   void toParseDown(TopSegmentIteratorConstPtr ts) const;
   BottomSegment* getBottomSegment();
   const BottomSegment* getBottomSegment() const;
   bool equals(BottomSegmentIteratorConstPtr other) const;
   bool hasChild(hal_size_t child) const;
   bool hasParseUp() const;
   bool inRange() const;

protected:
   HDF5BottomSegment _bottomSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable hal_bool_t _reversed;
};

inline bool HDF5BottomSegmentIterator::inRange() const
{
  return _bottomSegment._index >= 0 && _bottomSegment._index < 
     (hal_index_t)_bottomSegment._genome->_bottomArray.getSize();
}

}
#endif
