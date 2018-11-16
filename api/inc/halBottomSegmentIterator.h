/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENTITERATOR_H
#define _HALBOTTOMSEGMENTITERATOR_H

#include <cassert>
#include "halDefs.h"
#include "halBottomSegment.h"
#include "halSegmentIterator.h"
#include "halGenome.h"


namespace hal {

/** 
 * Interface for bottom segment iterator exposes the bottom segment
 * interface and some new methods for jumping around the genome.  
 * Always hidden in smart pointers in the public interface. 
 */
class BottomSegmentIterator : public virtual SegmentIterator
{
public:
    /** constructor */
    BottomSegmentIterator(BottomSegment* bottomSegment, 
                          hal_size_t startOffset = 0, 
                          hal_size_t endOffset = 0,
                          bool reversed = false) :
        SegmentIterator(startOffset, endOffset, reversed),
        _bottomSegment(bottomSegment) {
    }
    
    /** destructor */
    ~BottomSegmentIterator() {
    }

   /** Return a new copy of the iterator */
    BottomSegmentIteratorPtr clone() const;

   /** Copy an input iterator.  More efficient than the above methods
    * as no new iterator needs to be allocated 
    * @param botSegIt Iterator to copy */
    void copy(const BottomSegmentIteratorPtr& botSegIt);

   /** Move the iterator to the parent segment of a given iterator
    * @param topSegIt Iterator whose parent to move to */
    void toParent(const TopSegmentIteratorPtr& topSegIt); 

   /** Move the iterator down to the bottom segment containing the
    * start position of the given iterator in the same genome
    * @param topSegIt Top iterator to parse down on */
    void toParseDown(const TopSegmentIteratorPtr& topSegIt);

   /** DEPRECATED: FIXME: or not  */
    BottomSegment* getBottomSegment() {
        return  _bottomSegment.get();
    }

   /** DEPRECATED FIXME: or not */
    const BottomSegment* getBottomSegment() const {
        return  _bottomSegment.get();
    }

    /** return a pointer to the current BottomSegment (terse) */
    BottomSegment* bs() {
       return _bottomSegment.get();
     }

   /** return a pointer to the current BottomSegment (terse) */
   const BottomSegment* bs() const {
       return _bottomSegment.get();
   }

   /** Test equality with other iterator (current implementation does not
    * take into account reverse state or offsets -- FIXME: too review). 
    * FIXME merge with operator==?? 
    * @param other Iterator to test equality to */
    bool equals(const BottomSegmentIteratorPtr& other) const {
        assert(_bottomSegment->getGenome() == other->getGenome());
        return getArrayIndex() == other->getArrayIndex();
    }

    /* equality operator */
    bool operator==(const BottomSegmentIterator& other) const {
        assert(_bottomSegment->getGenome() == other.getBottomSegment()->getGenome());
        return getArrayIndex() == other.getArrayIndex();
    }

    /* inequality operator */
    bool operator!=(const BottomSegmentIterator& other) const {
        return !(*this == other);
    }
    
    // FIXME: document or change way getting segment works
    virtual Segment* getSegment() {
        return _bottomSegment.get();
    }
    virtual const Segment* getSegment() const {
        return _bottomSegment.get();
    }
    
    virtual void print(std::ostream& os) const;

private:
    virtual hal_size_t getNumSegmentsInGenome() const {
        return getGenome()->getNumBottomSegments();
    }
   BottomSegmentPtr _bottomSegment;
};

inline bool operator==(BottomSegmentIteratorPtr p1,
                       BottomSegmentIteratorPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(BottomSegmentIteratorPtr p1,
                       BottomSegmentIteratorPtr p2)
{
  return !(p1 == p2);
}
}


#endif
