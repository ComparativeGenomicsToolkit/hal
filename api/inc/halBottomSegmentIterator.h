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
class BottomSegmentIterator : public virtual BottomSegment,
                              public virtual SegmentIterator
{
public:
   /** Return a new copy of the iterator */
   virtual BottomSegmentIteratorPtr copy() = 0;

   /** Return a new copy of the iterator */
   virtual BottomSegmentIteratorConstPtr copy() const = 0;

   /** Copy an input iterator.  More efficient than the above methods
    * as no new iterator needs to be allocated 
    * @param ts Iterator to copy */
   virtual void copy(BottomSegmentIteratorConstPtr bs) const = 0;

   /** Move the iterator to the parent segment of a given iterator
    * @param ts Iterator whose parent to move to */
   virtual void toParent(TopSegmentIteratorConstPtr ts) const = 0; 

   /** Move the iterator down to the bottom segment containg the
    * start position of the given iterator in the same genome
    * @param ts Top iterator to parse down on */
   virtual void toParseDown(TopSegmentIteratorConstPtr ts) const = 0;

   /** DEPRECATED */
   virtual BottomSegment* getBottomSegment() = 0;

   /** DEPRECATED */
   virtual const BottomSegment* getBottomSegment() const = 0;

   /** Test equality with other iterator (current implementation does not
    * take into account reverse state or offsets -- too review)
    * @param other Iterator to test equality to */
   virtual bool equals(BottomSegmentIteratorConstPtr other) const = 0;

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
