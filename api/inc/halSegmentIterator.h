/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"

namespace hal {

/** 
 * Interface for general segment iterator
 */
class SegmentIterator
{
public:
   virtual ~SegmentIterator() = 0;
   virtual SegmentIteratorConstPtr parse() const = 0;
   virtual SegmentIteratorPtr parse() = 0;
   virtual std::vector<SegmentIteratorConstPtr> parseList() const = 0;
   virtual std::vector<SegmentIteratorPtr> parseList() = 0;
   virtual SegmentIteratorConstPtr parent() const = 0;
   virtual SegmentIteratorPtr parent() = 0;
   virtual SegmentIteratorConstPtr next() const = 0;
   virtual SegmentIteratorPtr next() = 0;
   virtual SegmentIteratorConstPtr prev() const = 0;
   virtual SegmentIteratorPtr prev() = 0;
   virtual SegmentIteratorConstPtr child(hal_offset_t) const = 0;
   virtual SegmentIteratorPtr child(hal_offset_t) = 0;
   virtual hal_offset_t leftOffset() const = 0;
   virtual hal_offset_t rightOffset() const = 0;
   virtual hal_bool_t amTop() const = 0;
   virtual hal_bool_t amBottom() const = 0;
   virtual SegmentConstPtr segment() const = 0;
   virtual SegmentPtr segment() = 0;
   virtual TopSegmentConstPtr topSegment() const = 0;
   virtual TopSegmentPtr topSegment() = 0;
   virtual BottomSegmentConstPtr bottomSegment() const = 0;
   virtual BottomSegmentPtr bototmSegment() = 0;
};



}
#endif
