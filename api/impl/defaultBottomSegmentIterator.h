/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _DEFAULTBOTTOMSEGMENTITERATOR_H
#define _DEFAULTBOTTOMSEGMENTITERATOR_H

#include "halBottomSegmentIterator.h"

namespace hal {

class DefaultBottomSegmentIterator : public BottomSegmentIterator
{
public:
   
   DefaultBottomSegmentIterator(BottomSegment* segment,
                             hal_size_t startOffset = 0, 
                             hal_size_t endOffset = 0,
                             bool inverted = false);
   virtual ~DefaultBottomSegmentIterator();
   
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
   
   // BOTTOM SEGMENT INTERFACE
   virtual hal_size_t getNumChildren() const;
   virtual hal_index_t getChildIndex(hal_size_t i) const;
   virtual hal_index_t getChildIndexG(const Genome* childGenome) const;
   virtual bool hasChild(hal_size_t child) const;
   virtual bool hasChildG(const Genome* childGenome) const;
   virtual void setChildIndex(hal_size_t i, hal_index_t childIndex);
   virtual bool getChildReversed(hal_size_t i) const;
   virtual void setChildReversed(hal_size_t child, bool isReversed);
   virtual hal_index_t getTopParseIndex() const;
   virtual void setTopParseIndex(hal_index_t parseIndex);
   virtual hal_offset_t getTopParseOffset() const;
   virtual bool hasParseUp() const;
   virtual hal_index_t getLeftChildIndex(hal_size_t i) const;
   virtual hal_index_t getRightChildIndex(hal_size_t i) const;

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

   // BOTTOM SEGMENT ITERATOR INTERFACE
   virtual BottomSegmentIteratorPtr copy();
   virtual BottomSegmentIteratorConstPtr copy() const;
   virtual void copy(BottomSegmentIteratorConstPtr bs) const;
   virtual void toParent(TopSegmentIteratorConstPtr ts) const; 
   virtual void toParseDown(TopSegmentIteratorConstPtr ts) const;
   virtual BottomSegment* getBottomSegment();
   virtual const BottomSegment* getBottomSegment() const;
   virtual bool equals(BottomSegmentIteratorConstPtr other) const;

protected:
   bool inRange() const;

protected:
   BottomSegmentPtr _bottomSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable bool _reversed;
};

inline bool DefaultBottomSegmentIterator::inRange() const
{
  return getArrayIndex() >= 0 && getArrayIndex() < 
     (hal_index_t)getGenome()->getNumBottomSegments();
}

}
#endif
