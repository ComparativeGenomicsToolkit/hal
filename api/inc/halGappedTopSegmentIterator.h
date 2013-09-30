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
 * Interface for Gepped Top Segment iterator.  Only used internally
 * and probably shouldn't be in public interface.
 */
class GappedTopSegmentIterator : virtual public GappedSegmentIterator
{
public:
   
   /** Return a copy of the iterator */
   virtual GappedTopSegmentIteratorPtr copy() = 0;

   /** Return a copy of the iterator */
   virtual GappedTopSegmentIteratorConstPtr copy() const = 0;

   /** Copy another iterator into the current iterator (more efficient
    * than above methods since no new iterators are created */
   virtual void copy(GappedTopSegmentIteratorConstPtr ts) const = 0;

   /** Move to child of given iterator 
    * @param bs Move to child of this iterator */
   virtual void toChild(GappedBottomSegmentIteratorConstPtr bs) const = 0;

   /** Test equality with other iterator 
    * @param other */
   virtual bool equals(GappedTopSegmentIteratorConstPtr other) const = 0;

   /** Test if iterator abuts other iterator */
   virtual bool adjacentTo(GappedTopSegmentIteratorConstPtr other) const = 0;

   /** Test if iterator has a parent */
   virtual bool hasParent() const = 0;

   /** Test if iterator has a duplicate */
   virtual bool hasNextParalogy() const = 0;

   /** Move to next paralogy */
   virtual void toNextParalogy() const = 0;

   /** Test if in reverse orientationt with respect to parent */
   virtual bool getParentReversed() const = 0;

   /** Return the leftmost segment of the iterator
    * (note that moving the returned iterator will corrupt the 
    * current gapped iterator.  this is a bug) */
   virtual TopSegmentIteratorConstPtr getLeft() const = 0;

   /** Return the rightmost segment of the iterator
    * (note that moving the returned iterator will corrupt the 
    * current gapped iterator.  this is a bug) */
   virtual TopSegmentIteratorConstPtr getRight() const = 0;

   /** Reset the gapped iterator.
    * @param ts This will be the left segment of the current iterator. The 
    * right segment will be extended as far as possible */
   virtual void setLeft(TopSegmentIteratorConstPtr ts) const = 0;

   /** For every set of paralogous top segments in a given genome, we identify a
    * unique segment as the canonical reference.  This is the segment that 
    * will be traversed when disabling duplication edges when mapping via 
    * the column iterator or mappedSegment interface.  It is also the segment
    * that is connected from its parent's down edge.*/
   virtual bool isCanonicalParalog() const = 0;

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

#endif
