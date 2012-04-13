/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENTITERATOR_H
#define _HALBOTTOMSEGMENTITERATOR_H

#include "halDefs.h"
#include "halBottomSegment.h"
#include "halSegmentIterator.h"

namespace hal {

/** 
 * Interface for bottom segment iterator exposes the bottom segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class BottomSegmentIterator : public SegmentIterator
{
public:
   virtual BottomSegmentIteratorPtr copy() = 0;
   virtual BottomSegmentIteratorConstPtr copy() const = 0;
   virtual void toParent(TopSegmentIteratorConstPtr ts) const = 0; 
   virtual void toParseDown(TopSegmentIteratorConstPtr bs) const = 0;
   virtual BottomSegment* getBottomSegment() = 0;
   virtual const BottomSegment* getBottomSegment() const = 0;

protected:
   friend class counted_ptr<BottomSegmentIterator>;
   friend class counted_ptr<const BottomSegmentIterator>;
   virtual ~BottomSegmentIterator() = 0;
};

inline BottomSegmentIterator::~BottomSegmentIterator() {}

}
#endif
