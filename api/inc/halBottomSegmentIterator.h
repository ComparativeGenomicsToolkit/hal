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
class BottomSegmentIterator : public virtual BottomSegment,
                              public virtual SegmentIterator
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
    void copy(BottomSegmentIteratorPtr botSegIt);

   /** Move the iterator to the parent segment of a given iterator
    * @param topSegIt Iterator whose parent to move to */
    void toParent(TopSegmentIteratorPtr topSegIt); 

   /** Move the iterator down to the bottom segment containing the
    * start position of the given iterator in the same genome
    * @param topSegIt Top iterator to parse down on */
    void toParseDown(TopSegmentIteratorPtr topSegIt);

   /** DEPRECATED */
    BottomSegment* getBottomSegment() {
        return  _bottomSegment.get();
    }

   /** DEPRECATED */
    const BottomSegment* getBottomSegment() const {
        return  _bottomSegment.get();
    }

   /** Test equality with other iterator (current implementation does not
    * take into account reverse state or offsets -- FIXME: too review). 
    * FIXME merge with operator==?? 
    * @param other Iterator to test equality to */
    bool equals(BottomSegmentIteratorPtr other) const {
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
    
    // SEGMENT INTERFACE OVERRIDE
    virtual void print(std::ostream& os) const;

   // BOTTOM SEGMENT INTERFACE
    virtual hal_size_t getNumChildren() const {
        return _bottomSegment->getNumChildren();
    }
    virtual hal_index_t getChildIndex(hal_size_t i) const {
          return _bottomSegment->getChildIndex(i);
    }
    virtual hal_index_t getChildIndexG(const Genome* childGenome) const {
          return _bottomSegment->getChildIndexG(childGenome);
    }
    virtual bool hasChild(hal_size_t child) const {
        return _bottomSegment->hasChild(child);
    }
    virtual bool hasChildG(const Genome* childGenome) const {
        return _bottomSegment->hasChildG(childGenome);

    }
    virtual void setChildIndex(hal_size_t i, hal_index_t childIndex) {
        _bottomSegment->setChildIndex(i, childIndex);
    }
    virtual bool getChildReversed(hal_size_t i) const {
          return _bottomSegment->getChildReversed(i);
    }
    virtual void setChildReversed(hal_size_t child, bool isReversed) {
        _bottomSegment->setChildReversed(child, isReversed);
    }
    virtual hal_index_t getTopParseIndex() const {
        return _bottomSegment->getTopParseIndex();
    }
    virtual void setTopParseIndex(hal_index_t parseIndex) {
        _bottomSegment->setTopParseIndex(parseIndex);
    }
    virtual hal_offset_t getTopParseOffset() const {
        return _bottomSegment->getTopParseOffset();
    }
    virtual bool hasParseUp() const {
        return _bottomSegment->hasParseUp();
    }
    virtual hal_index_t getLeftChildIndex(hal_size_t i) const {
          assert(_startOffset == 0 && _endOffset == 0);
          return _bottomSegment->getLeftChildIndex(i);
    }
    virtual hal_index_t getRightChildIndex(hal_size_t i) const {
        assert(_startOffset == 0 && _endOffset == 0);
        return _bottomSegment->getRightChildIndex(i);
    }
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
