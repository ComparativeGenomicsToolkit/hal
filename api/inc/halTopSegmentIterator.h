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
#include "halGenome.h"
#include "halSegmentIterator.h"

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
                      bool inverted = false);

    /* destructor */
    virtual ~TopSegmentIterator() {
    }
   
   /** Return a new copy of the iterator */
    TopSegmentIteratorPtr copy();

   /** Return a new copy of the iterator */
     TopSegmentIteratorConstPtr copy() const;

   /** Copy an input iterator.  More efficient than the above methods
    * as no new iterator needs to be allocated 
    * @param ts Iterator to copy */
    void copy(TopSegmentIteratorConstPtr ts) const;

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param child Index of child in bottom segment's genome */
    void toChild(BottomSegmentIteratorConstPtr bs, 
                 hal_size_t child) const;

   /** Move the iterator to the child of a given bottom segment
    * @param bs Bottom segment whose child will be moved to
    * @param childGenome genome of child in bottom segment */
    void toChildG(BottomSegmentIteratorConstPtr bs, 
                  const Genome* childGenome) const;
   
   /** Given a bottom segment, move to the top segment that contains
    * its start position.  The genome remains unchanged.  The iterator
    * will be sliced accordingly (reversed state also taken into account)
    * @param bs Bottom segment to parse up from */
    void toParseUp(BottomSegmentIteratorConstPtr bs) const;

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
    * @param other Iterator to test equality to */
    bool equals(TopSegmentIteratorConstPtr other) const;

   /** Move iterator to next paralgous segment.  Iterator will be reversed
   * if the next segment is in a different orientation wrt their common
   * parent */
    void toNextParalogy() const;

    // FIXME: document or change way getting segment works
    virtual SegmentPtr getSegment();
    virtual SegmentConstPtr getSegment() const;
    
    // SEGMENT INTERFACE OVERRIDE
    virtual void print(std::ostream& os) const;
   // TOP SEGMENT INTERFACE
   virtual hal_index_t getParentIndex() const;
   virtual bool hasParent() const;
   virtual void setParentIndex(hal_index_t parIdx);
   virtual bool getParentReversed() const;
   virtual void setParentReversed(bool isReversed);
   virtual hal_index_t getBottomParseIndex() const;
   virtual void setBottomParseIndex(hal_index_t botParseIdx);
   virtual hal_offset_t getBottomParseOffset() const;
   virtual bool hasParseDown() const;
   virtual hal_index_t getNextParalogyIndex() const;
   virtual bool hasNextParalogy() const;
   virtual void setNextParalogyIndex(hal_index_t parIdx);
   virtual hal_index_t getLeftParentIndex() const;
   virtual hal_index_t getRightParentIndex() const;
   virtual bool isCanonicalParalog() const;


private:
   virtual hal_size_t getNumSegmentsInGenome() const;

   TopSegmentPtr _topSegment;

   friend class counted_ptr<TopSegmentIterator>;
   friend class counted_ptr<const TopSegmentIterator>;
};

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

inline TopSegmentIterator::TopSegmentIterator(TopSegment* topSegment, 
                                             hal_offset_t startOffset, 
                                             hal_offset_t endOffset,
                                             bool reversed) :
    SegmentIterator(startOffset, endOffset, reversed),
    _topSegment(topSegment)
{

}

inline SegmentPtr TopSegmentIterator::getSegment()
{
  return _topSegment;
}

inline SegmentConstPtr TopSegmentIterator::getSegment() const
{
  return _topSegment;
}

inline hal_size_t TopSegmentIterator::getNumSegmentsInGenome() const
{
  return getGenome()->getNumTopSegments();
}


}


#endif
