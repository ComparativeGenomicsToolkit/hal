/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDSEGMENTITERATOR_H
#define _HALGAPPEDSEGMENTITERATOR_H

#include "halSegmentIterator.h"

namespace hal {

/** 
 * Interface for general gappedSegment iterator.  Behaves like
 * a regular iterator, but operates on a linear sequence of 
 * segments that are consistent modulo gaps.  
 */
class GappedSegmentIterator : public SegmentIterator
{
public:
   
   virtual hal_size_t getGapThreshold() const = 0;

   virtual hal_size_t getChildIndex() const = 0;
   
   virtual hal_size_t getNumSegments() const = 0;

   virtual hal_size_t getNumGaps() const = 0;
   
   virtual hal_size_t getNumGapBases() const = 0;
   
   virtual bool isLast() const = 0;
   
   virtual bool isFirst() const = 0;

   virtual hal_index_t getLeftArrayIndex() const = 0;

   virtual hal_index_t getRightArrayIndex() const = 0;

   virtual const Sequence* getSequence() const = 0;

protected:
   friend class counted_ptr<GappedSegmentIterator>;
   friend class counted_ptr<const GappedSegmentIterator>;
   virtual ~GappedSegmentIterator() = 0;
};

inline GappedSegmentIterator::~GappedSegmentIterator() {}

}
#endif
