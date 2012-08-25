/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDEFAULTGAPPEDTOPSEGMENTITERATOR_H
#define _HALDEFAULTGAPPEDTOPSEGMENTITERATOR_H

#include "halGappedTopSegmentIterator.h"

namespace hal {

/** 
 * Default (implementation independent) implementation of 
 * GappedTopSegmentIterator interface 
 */
class DefaultGappedTopSegmentIterator : public GappedTopSegmentIterator
{
public:

   DefaultGappedTopSegmentIterator(TopSegmentIteratorConstPtr left,
                                   hal_size_t gapThreshold);

   ~DefaultGappedTopSegmentIterator();

   // Gpped Segment Iterator methods
   hal_size_t getGapThreshold() const;
   hal_size_t getChildIndex() const;
   hal_size_t getNumSegments() const;
   hal_size_t getNumGaps() const;
   hal_size_t getNumGapBases() const;
   bool isLast() const;
   bool isFirst() const;
   hal_index_t getLeftArrayIndex() const;
   hal_index_t getRightArrayIndex() const;
   const Sequence* getSequence() const;
   
   // Segment Iterator methods
   void toLeft(hal_index_t leftCutoff) const;
   void toRight(hal_index_t rightCutoff) const;
   void toReverse() const;
   void toSite(hal_index_t position, bool slice = true) const;
   bool hasNextParalogy() const;
   void toNextParalogy() const;
   hal_offset_t getStartOffset() const;
   hal_offset_t getEndOffset() const;
   void slice(hal_offset_t startOffset,
              hal_offset_t endOffset) const;
   hal_index_t getStartPosition() const;
   hal_size_t getLength() const;
   hal_bool_t getReversed() const;
   void getString(std::string& outString) const;
   bool leftOf(hal_index_t genomePos) const;
   bool rightOf(hal_index_t genomePos) const;
   bool overlaps(hal_index_t genomePos) const;

   // GappedTopSegmentIterator Interface Methods
   GappedTopSegmentIteratorPtr copy();
   GappedTopSegmentIteratorConstPtr copy() const;
   void copy(GappedTopSegmentIteratorConstPtr ts) const;
   void toChild(GappedBottomSegmentIteratorConstPtr bs) const;
   bool equals(GappedTopSegmentIteratorConstPtr other) const;
   bool adjacentTo(GappedTopSegmentIteratorConstPtr other) const;
   bool hasParent() const;
   bool getParentReversed() const;
   TopSegmentIteratorConstPtr getLeft() const;
   TopSegmentIteratorConstPtr getRight() const;
   void setLeft(TopSegmentIteratorConstPtr ts) const;


private:
   
   bool compatible(TopSegmentIteratorConstPtr left,
                   TopSegmentIteratorConstPtr right) const;

   void extendRight() const;
   void extendLeft() const;

   void toLeftNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toRightNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toLeftNextUngapped(TopSegmentIteratorConstPtr ts) const;
   void toRightNextUngapped(TopSegmentIteratorConstPtr ts) const;
   
   // keep convention of other iterators where const-ness only applies
   // to the database and not the iterator...
   mutable TopSegmentIteratorConstPtr _left;
   mutable TopSegmentIteratorConstPtr _right;
   mutable BottomSegmentIteratorConstPtr _leftParent;
   mutable BottomSegmentIteratorConstPtr _rightParent;
   mutable TopSegmentIteratorConstPtr _leftDup;
   mutable TopSegmentIteratorConstPtr _rightDup;
   mutable TopSegmentIteratorConstPtr _temp;
   mutable TopSegmentIteratorConstPtr _temp2;
   mutable hal_size_t _childIndex;
   mutable hal_size_t _gapThreshold;            
  
};


}
#endif
