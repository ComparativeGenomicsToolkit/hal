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
   void toLeft() const;
   void toRight() const;
   void toNextParalogy() const;
   hal_offset_t getStartOffset() const;
   hal_offset_t getEndOffset() const;
   hal_size_t getLength() const;
   hal_bool_t getReversed() const;
   void getSequence(std::string& outSequence);

   //TOP ITERATOR METHODS
   TopSegmentIteratorPtr copy();
   TopSegmentIteratorConstPtr copy() const;
   void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const;
   void toParseUp(BottomSegmentIteratorConstPtr bs) const;
   TopSegment* getTopSegment();
   const TopSegment* getTopSegment() const;

protected:
   HDF5TopSegment _topSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable hal_bool_t _reversed;
};

}
#endif
