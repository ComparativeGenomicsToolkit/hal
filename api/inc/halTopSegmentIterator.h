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
#include "halSegmentIterator.h"
#include "halGenome.h"

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
    /* constructor */
   TopSegmentIterator(TopSegment* topSegment,
                      hal_offset_t startOffset = 0, 
                      hal_offset_t endOffset = 0,
                      bool reversed = false) :
       SegmentIterator(startOffset, endOffset, reversed),
       _topSegment(topSegment) {
    
   }
      

    /* destructor */
    virtual ~TopSegmentIterator() {
    }
   
   /** Return a new copy of the iterator */
     TopSegmentIteratorPtr clone() const;

   /** Copy an input iterator.  More efficient than the above methods
    * as no new iterator needs to be allocated 
    * @param ts Iterator to copy */
    void copy(TopSegmentIteratorPtr ts);

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param child Index of child in bottom segment's genome */
    void toChild(BottomSegmentIteratorPtr bs, 
                 hal_size_t child);

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param childGenome genome of child in bottom segment */
    void toChildG(BottomSegmentIteratorPtr bs, 
                  const Genome* childGenome);
   
   /** Given a bottom segment, move to the top segment that contains
    * its start position.  The genome remains unchanged.  The iterator
    * will be sliced accordingly (reversed state also taken into account)
    * @param bs Bottom segment to parse up from */
    void toParseUp(BottomSegmentIteratorPtr bs);

   /** DEPRECATED */
    TopSegment* getTopSegment() {
        return _topSegment.get();
    }

   /** DEPRECATED */
    const TopSegment* getTopSegment() const  {
        return _topSegment.get();
   }

   /** Test equality with other iterator (current implementation does not
    * take into account reverse state or offsets -- too review)
    * FIXME merge with operator==?? 
    * @param other Iterator to test equality to */
    bool equals(TopSegmentIteratorPtr other) const {
        assert(getGenome() == other->getGenome());
        return getArrayIndex() == other->getArrayIndex();
    }


    /* equality operator */
    bool operator==(const TopSegmentIterator& other) const {
        assert(_topSegment->getGenome() == other.getTopSegment()->getGenome());
        return getArrayIndex() == other.getArrayIndex();
    }

    /* inequality operator */
    bool operator!=(const TopSegmentIterator& other) const {
        return !(*this == other);
    }
    
   /** Move iterator to next paralgous segment.  Iterator will be reversed
   * if the next segment is in a different orientation wrt their common
   * parent */
    void toNextParalogy();

    // FIXME: document or change way getting segment works
    virtual Segment* getSegment() {
        return _topSegment.get();
    }
    virtual const Segment* getSegment() const {
        return _topSegment.get();
    }
    
    // SEGMENT INTERFACE OVERRIDE
    virtual void print(std::ostream& os) const;
   // TOP SEGMENT INTERFACE
    virtual hal_index_t getParentIndex() const {
  return _topSegment->getParentIndex();
    }
    virtual bool hasParent() const {
        assert(inRange());
        return _topSegment->getParentIndex() != NULL_INDEX;
    }
    virtual void setParentIndex(hal_index_t parIdx) {
        _topSegment->setParentIndex(parIdx);
    }
    virtual bool getParentReversed() const {
        return _topSegment->getParentReversed();
    }
    virtual void setParentReversed(bool isReversed) {
        _topSegment->setParentReversed(isReversed);
    }
    virtual hal_index_t getBottomParseIndex() const {
        return _topSegment->getBottomParseIndex();
    }
    virtual void setBottomParseIndex(hal_index_t botParseIdx) {
        _topSegment->setBottomParseIndex(botParseIdx);
    }
    virtual hal_offset_t getBottomParseOffset() const {
        return _topSegment->getBottomParseOffset();
    }
    virtual bool hasParseDown() const {
        assert (inRange() == true);
        assert (_topSegment->getBottomParseIndex() == NULL_INDEX ||
                _topSegment->getGenome()->getNumChildren() > 0);
        return _topSegment->getBottomParseIndex() != NULL_INDEX;
    }
    virtual hal_index_t getNextParalogyIndex() const {
        return _topSegment->getNextParalogyIndex();
    }
    virtual bool hasNextParalogy() const {
        return _topSegment->hasNextParalogy();
    }
    virtual void setNextParalogyIndex(hal_index_t parIdx) {
        _topSegment->setNextParalogyIndex(parIdx);
    }
    virtual hal_index_t getLeftParentIndex() const {
        return _topSegment->getLeftParentIndex();
    }
    virtual hal_index_t getRightParentIndex() const {
        return _topSegment->getRightParentIndex();
    }
    virtual bool isCanonicalParalog() const {
        return _topSegment->isCanonicalParalog();
    }


private:
    hal_size_t getNumSegmentsInGenome() const {
        return getGenome()->getNumTopSegments();
    }

   TopSegmentPtr _topSegment;
};

inline bool operator==(TopSegmentIteratorPtr p1,
                       TopSegmentIteratorPtr p2) 
{
  if (p1.get() == NULL || p2.get() == NULL)
  {
    return p1.get() == NULL && p2.get() == NULL;
  }
  return p1->equals(p2);
}

inline bool operator!=(TopSegmentIteratorPtr p1,
                       TopSegmentIteratorPtr p2)
{
  return !(p1 == p2);
}
}

#endif
