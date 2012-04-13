/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for general segment iterator.  Common functionality
 * for top and bottom iterators.  A segment iterator implements 
 * the segment as well as the iterator interface. 
 */
class SegmentIterator : public Segment
{
public:

   virtual SegmentIteratorPtr copy() = 0;
   virtual SegmentIteratorConstPtr copy() const = 0;
   virtual void toLeft() const = 0;
   virtual void toRight() const = 0;
   virtual hal_offset_t getStartOffset() const = 0;
   virtual hal_offset_t getEndOffset() const = 0;

protected:
   friend class counted_ptr<SegmentIterator>;
   friend class counted_ptr<const SegmentIterator>;
   virtual ~SegmentIterator() = 0;
};

inline SegmentIterator::~SegmentIterator() {}

}
#endif
