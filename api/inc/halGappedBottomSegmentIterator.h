/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDBOTTOMSEGMENTITERATOR_H
#define _HALGAPPEDBOTTOMSEGMENTITERATOR_H

#include <iostream>
#include "halDefs.h"
#include "halBottomSegmentIterator.h"
#include "halGappedSegmentIterator.h"

namespace hal {

/** 
 * Interface for gappedBottom segment iterator exposes the gappedBottom segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class GappedBottomSegmentIterator : public GappedSegmentIterator
{
public:
   
   virtual GappedBottomSegmentIteratorPtr copy() = 0;
   virtual GappedBottomSegmentIteratorConstPtr copy() const = 0;
   virtual void copy(GappedBottomSegmentIteratorConstPtr ts) const = 0;

   virtual void toParent(GappedTopSegmentIteratorConstPtr ts) const = 0;
   virtual bool equals(GappedBottomSegmentIteratorConstPtr other) const = 0;
   virtual bool adjacentTo(GappedBottomSegmentIteratorConstPtr other) const= 0;
   virtual bool hasChild() const = 0;
   virtual bool getChildReversed() const = 0;

   virtual BottomSegmentIteratorConstPtr getLeft() const = 0;
   virtual BottomSegmentIteratorConstPtr getRight() const = 0;

protected:
   friend class counted_ptr<GappedBottomSegmentIterator>;
   friend class counted_ptr<const GappedBottomSegmentIterator>;
   virtual ~GappedBottomSegmentIterator() = 0;
};

inline GappedBottomSegmentIterator::~GappedBottomSegmentIterator() {}

inline bool operator==(GappedBottomSegmentIteratorConstPtr p1,
                       GappedBottomSegmentIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(GappedBottomSegmentIteratorConstPtr p1,
                       GappedBottomSegmentIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}


#ifndef NDEBUG
#include "halGenome.h"
namespace hal {
inline std::ostream& operator<<(std::ostream& os, 
                                GappedBottomSegmentIteratorConstPtr gbs)
{
  os << "Th=" << gbs->getGapThreshold() << " chIdx=" << gbs->getChildIndex() 
     << "\nLEFT: " << gbs->getLeft() << "\nRIGHT: " << gbs->getRight();
  return os;
}
}
#endif
#endif
