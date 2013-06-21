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
 * Interface for Gepped Bottom Segment iterator.  Only used internally
 * and probably shouldn't be in public interface.
 */
class GappedBottomSegmentIterator : virtual public GappedSegmentIterator
{
public:
   
   /** Return a copy of the iterator */
   virtual GappedBottomSegmentIteratorPtr copy() = 0;

   /** Return a copy of the iterator */
   virtual GappedBottomSegmentIteratorConstPtr copy() const = 0;

   /** Copy another iterator into the current iterator (more efficient
    * than above methods since no new iterators are created */
   virtual void copy(GappedBottomSegmentIteratorConstPtr ts) const = 0;

   /** Move to parent */
   virtual void toParent(GappedTopSegmentIteratorConstPtr ts) const = 0;

   /** Test equality with other iterator 
    * @param other */
   virtual bool equals(GappedBottomSegmentIteratorConstPtr other) const = 0;

   /** Test if iterator abuts other iterator */
   virtual bool adjacentTo(GappedBottomSegmentIteratorConstPtr other) const= 0;

   /** Test if iterator has a child */
   virtual bool hasChild() const = 0;

   /** Test if in reverse orientationt with respect to child */
   virtual bool getChildReversed() const = 0;

   /** Return the rightmost segment of the iterator
    * (note that moving the returned iterator will corrupt the 
    * current gapped iterator.  this is a bug) */
   virtual BottomSegmentIteratorConstPtr getLeft() const = 0;

   /** Reset the gapped iterator.
    * @param ts This will be the left segment of the current iterator. The 
    * right segment will be extended as far as possible */
   virtual BottomSegmentIteratorConstPtr getRight() const = 0;

   /** Reset the gapped iterator.
    * @param ts This will be the left segment of the current iterator. The 
    * right segment will be extended as far as possible */
   virtual void setLeft(BottomSegmentIteratorConstPtr bs) const = 0;


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


#endif
