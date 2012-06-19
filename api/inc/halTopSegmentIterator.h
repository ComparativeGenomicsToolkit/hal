/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENTITERATOR_H
#define _HALTOPSEGMENTITERATOR_H

#include <iostream>
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
   virtual void toSite(hal_index_t position, bool slice = true) const = 0;
   virtual TopSegment* getTopSegment() = 0;
   virtual const TopSegment* getTopSegment() const = 0;
   virtual bool equals(TopSegmentIteratorConstPtr other) const = 0;
   virtual bool hasParent() const = 0;
   virtual bool hasParseDown() const = 0;

protected:
   friend class counted_ptr<TopSegmentIterator>;
   friend class counted_ptr<const TopSegmentIterator>;
   virtual ~TopSegmentIterator() = 0;
};

inline TopSegmentIterator::~TopSegmentIterator() {}

inline bool operator==(TopSegmentIteratorConstPtr p1,
                       TopSegmentIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(TopSegmentIteratorConstPtr p1,
                       TopSegmentIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}


#ifndef NDEBUG
#include "halGenome.h"
namespace hal {
inline std::ostream& operator<<(std::ostream& os, 
                                TopSegmentIteratorConstPtr ts)
{
  const TopSegment* tseg = ts->getTopSegment();
  const Genome* genome = tseg->getGenome();
  os << "topIterator: ";
  os << "Genome=" << genome->getName();
  os << " idx=" << tseg->getArrayIndex();
  if (tseg->getArrayIndex() >= 0 && 
      tseg->getArrayIndex() < (hal_index_t)genome->getNumTopSegments())
  {
    os << " segLen=" << tseg->getLength();
    os << " segStart=" << tseg->getStartPosition() << std::endl;

    os << " offset=[" << ts->getStartOffset() << "," 
       << ts->getEndOffset() << "]";
    os << " start=" << ts->getStartPosition();
    os << " len=" << ts->getLength();
    os << " rev=" << ts->getReversed();
  }
  return os;
}
}
#endif

#endif
