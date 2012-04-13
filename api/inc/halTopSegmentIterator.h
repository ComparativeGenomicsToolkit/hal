/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENTITERATOR_H
#define _HALTOPSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSegmentIterator.h"
#include "halTopSegment.h"

namespace hal {

/** 
 * Interface for top segment iterator exposes the top segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class TopSegmentIterator : public SegmentIterator
{
public:
   virtual TopSegmentIteratorPtr copy() = 0;
   virtual TopSegmentIteratorConstPtr copy() const = 0;
   virtual void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const = 0;
   virtual void toParseUp(BottomSegmentIteratorConstPtr bs) const = 0;
   virtual TopSegment* getTopSegment() = 0;
   virtual const TopSegment* getTopSegment() const = 0;

protected:
   friend class counted_ptr<TopSegmentIterator>;
   friend class counted_ptr<const TopSegmentIterator>;
   virtual ~TopSegmentIterator() = 0;
};

inline TopSegmentIterator::~TopSegmentIterator() {}


}
#endif
