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
   virtual void copy(BottomSegmentIteratorConstPtr bs) const = 0;
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

#ifndef NDEBUG
#include "halGenome.h" 
namespace hal {
inline std::ostream& operator<<(std::ostream& os, 
                                BottomSegmentIteratorConstPtr bs)
{
  const BottomSegment* bseg = bs->getBottomSegment();
  const Genome* genome = bseg->getGenome();
  os << "bottomIterator: ";
  os << "Genome=" << genome->getName();
  os << " Seq=" << bseg->getSequence()->getName();
  os << " idx=" << bseg->getArrayIndex();
  if (bseg->getArrayIndex() >= 0 && 
      bseg->getArrayIndex() < (hal_index_t)genome->getNumBottomSegments())
  {
    os << " segLen=" << bseg->getLength();
    os << " segStart=" << bseg->getStartPosition() << std::endl;
    
    os << " offset=[" << bs->getStartOffset() << "," 
       << bs->getEndOffset() << "]";
    os << " start=" << bs->getStartPosition();
    os << " len=" << bs->getLength();
    os << " rev=" << bs->getReversed();
    for (hal_size_t i = 0; i < bseg->getNumChildren(); ++i)
    {
      os << " childIdx[" << i << "]=" << bseg->getChildIndex(i);
    }
  }
  return os;
}
}
#endif


#endif
