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
   void toLeft() const;
   void toRight() const;
   void toReverse() const;
   void toNextParalogy() const;
   hal_offset_t getStartOffset() const;
   hal_offset_t getEndOffset() const;
   hal_size_t getLength() const;
   hal_bool_t getReversed() const;
   void getString(std::string& outString) const;

   //BOTTOM ITERATOR METHODS
   BottomSegmentIteratorPtr copy();
   BottomSegmentIteratorConstPtr copy() const;
   void toParent(TopSegmentIteratorConstPtr ts) const; 
   void toParseDown(TopSegmentIteratorConstPtr ts) const;
   BottomSegment* getBottomSegment();
   const BottomSegment* getBottomSegment() const;
   bool equals(BottomSegmentIteratorConstPtr other) const;
   
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
