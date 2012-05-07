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
   virtual bool equals(BottomSegmentIteratorConstPtr other) const = 0;
   virtual bool hasChild(hal_size_t child) const = 0;
   virtual bool hasParseUp() const = 0;

protected:
   friend class counted_ptr<BottomSegmentIterator>;
   friend class counted_ptr<const BottomSegmentIterator>;
   virtual ~BottomSegmentIterator() = 0;
};

inline BottomSegmentIterator::~BottomSegmentIterator() {}

inline bool operator==(BottomSegmentIteratorConstPtr p1,
                       BottomSegmentIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(BottomSegmentIteratorConstPtr p1,
                       BottomSegmentIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}
#endif
