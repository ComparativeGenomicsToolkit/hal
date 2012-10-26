/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTTOPSEGMENTITERATOR_H
#define _DEFAULTTOPSEGMENTITERATOR_H

#include "halTopSegmentIterator.h"

namespace hal {


class DefaultTopSegmentIterator : public TopSegmentIterator
{
public:
   DefaultTopSegmentIterator(TopSegment* topSegment,
                          hal_offset_t startOffset = 0, 
                          hal_offset_t endOffset = 0,
                          bool inverted = false);
   virtual ~DefaultTopSegmentIterator();
   
   // SEGMENT INTERFACE
   virtual void setArrayIndex(Genome* genome, 
                              hal_index_t arrayIndex);
   virtual void setArrayIndex(const Genome* genome, 
                              hal_index_t arrayIndex) const;
   virtual const Genome* getGenome() const;
   virtual Genome* getGenome();
   virtual const Sequence* getSequence() const;
   virtual Sequence* getSequence();
   virtual hal_index_t getStartPosition() const;
   virtual hal_index_t getEndPosition() const;
   virtual hal_size_t getLength() const;
   virtual void getString(std::string& outString) const;
   virtual void setCoordinates(hal_index_t startPos, hal_size_t length);
   virtual hal_index_t getArrayIndex() const;
   virtual bool leftOf(hal_index_t genomePos) const;
   virtual bool rightOf(hal_index_t genomePos) const;
   virtual bool overlaps(hal_index_t genomePos) const;
   virtual bool isFirst() const;
   virtual bool isLast() const;
   virtual bool isMissingData(double nThreshold) const;

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

   // SEGMENT ITERATOR INTERFACE 
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   virtual void toReverse() const;
   virtual void toSite(hal_index_t position, bool slice = true) const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset ) const;
   virtual bool getReversed() const;

   // TOP SEGMENT ITERATOR INTERFACE 
   virtual TopSegmentIteratorPtr copy();
   virtual TopSegmentIteratorConstPtr copy() const;
   virtual void copy(TopSegmentIteratorConstPtr ts) const;
   virtual void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const;
   virtual void toChildG(BottomSegmentIteratorConstPtr bs, 
                         const Genome* childGenome) const;
   virtual void toParseUp(BottomSegmentIteratorConstPtr bs) const;
   virtual TopSegment* getTopSegment();
   virtual const TopSegment* getTopSegment() const;
   virtual bool equals(TopSegmentIteratorConstPtr other) const;
   virtual void toNextParalogy() const;

protected:
   bool inRange() const;
   
protected:
   TopSegmentPtr _topSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable bool _reversed;
};

inline bool DefaultTopSegmentIterator::inRange() const
{
  return getArrayIndex() >= 0 && getArrayIndex() < 
     (hal_index_t)getGenome()->getNumTopSegments();
}

}
#endif
