/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDTOPSEGMENTITERATOR_H
#define _HALGAPPEDTOPSEGMENTITERATOR_H

#include <iostream>
#include "halDefs.h"
#include "halTopSegmentIterator.h"
#include "halGappedSegmentIterator.h"

namespace hal {

/** 
 * Interface for gappedTop segment iterator exposes the gappedTop segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class GappedTopSegmentIterator : public GappedSegmentIterator
{
public:
   
   virtual GappedTopSegmentIteratorPtr copy() = 0;
   virtual GappedTopSegmentIteratorConstPtr copy() const = 0;
   virtual void copy(GappedTopSegmentIteratorConstPtr ts) const = 0;

   virtual void toChild(GappedBottomSegmentIteratorConstPtr bs) const = 0;
   virtual bool equals(GappedTopSegmentIteratorConstPtr other) const = 0;
   virtual bool hasParent() const = 0;

   virtual TopSegmentIteratorConstPtr getLeft() const = 0;
   virtual TopSegmentIteratorConstPtr getRight() const = 0;

protected:
   friend class counted_ptr<GappedTopSegmentIterator>;
   friend class counted_ptr<const GappedTopSegmentIterator>;
   virtual ~GappedTopSegmentIterator() = 0;
};

inline GappedTopSegmentIterator::~GappedTopSegmentIterator() {}

inline bool operator==(GappedTopSegmentIteratorConstPtr p1,
                       GappedTopSegmentIteratorConstPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(GappedTopSegmentIteratorConstPtr p1,
                       GappedTopSegmentIteratorConstPtr p2)
{
  return !(p1 == p2);
}

}


#ifndef NDEBUG
#include "halGenome.h"
namespace hal {
inline std::ostream& operator<<(std::ostream& os, 
                                GappedTopSegmentIteratorConstPtr gts)
{
  os << "first\n" << gts->getLeft() << "\nlast" << gts->getRight();
  return os;
}
}
#endif

#endif
