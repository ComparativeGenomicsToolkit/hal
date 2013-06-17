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
class TopSegmentIterator : public virtual TopSegment,
                           public virtual SegmentIterator
{
public:
   /** Return a new copy of the iterator */
   virtual TopSegmentIteratorPtr copy() = 0;

   /** Return a new copy of the iterator */
   virtual TopSegmentIteratorConstPtr copy() const = 0;

   /** Copy an input iterator.  More efficient than the above methods
    * as no new iterator needs to be allocated 
    * @param ts Iterator to copy */
   virtual void copy(TopSegmentIteratorConstPtr ts) const = 0;

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param child Index of child in bottom segment's genome */
   virtual void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const = 0;

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param childGenome genome of child in bottom segment */
   virtual void toChildG(BottomSegmentIteratorConstPtr bs, 
                         const Genome* childGenome) const = 0;
   
   /** Given a bottom segment, move to the top segment that contains
    * its start position.  The genome remains unchanged.  The iterator
    * will be sliced accordingly (reversed state also taken into account)
    * @param bs Bottom segment to parse up from */
   virtual void toParseUp(BottomSegmentIteratorConstPtr bs) const = 0;

   /** DEPRECATED */
   virtual TopSegment* getTopSegment() = 0;

   /** DEPRECATED */
   virtual const TopSegment* getTopSegment() const = 0;

   /** Test equality with other iterator (current implementation does not
    * take into account reverse state or offsets -- too review)
    * @param other Iterator to test equality to */
   virtual bool equals(TopSegmentIteratorConstPtr other) const = 0;

   /** Move iterator to next paralgous segment.  Iterator will be reversed
   * if the next segment is in a different orientation wrt their common
   * parent */
   virtual void toNextParalogy() const = 0;

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


#endif
