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

class TopSegmentIterator;
class TopSegmentConstIterator;

/** 
 * Interface for bottom segment iterator exposes the bottom segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class BottomSegmentIterator : public BottomSegment, public SegmentIterator
{
public:
   virtual void toParent(TopSegmentConstIterator ts) const = 0; 
   virtual void toParseDown(TopSegmentConstIterator bs) const = 0;

protected:
   friend class counted_ptr<BottomSegmentIterator>;
   friend class counted_ptr<const BottomSegmentIterator>;
   virtual ~BottomSegmentIterator() = 0;
};

inline BottomSegmentIterator::~BottomSegmentIterator() {}

}
#endif
