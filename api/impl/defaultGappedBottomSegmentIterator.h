/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDEFAULTGAPPEDBOTTOMSEGMENTITERATOR_H
#define _HALDEFAULTGAPPEDBOTTOMSEGMENTITERATOR_H

#include "halGappedBottomSegmentIterator.h"

namespace hal {

/** 
 * Default (implementation independent) implementation of 
 * GappedBottomSegmentIterator interface 
 */
class DefaultGappedBottomSegmentIterator : public GappedBottomSegmentIterator
{
public:

   DefaultGappedBottomSegmentIterator(BottomSegmentIteratorConstPtr left,
                                      hal_size_t childIndex,
                                      hal_size_t gapThreshold);

   ~DefaultGappedBottomSegmentIterator();

   // Gpped Segment Iterator methods
   hal_size_t getGapThreshold() const;
   hal_size_t getChildIndex() const;
   hal_size_t getNumSegments() const;
   hal_size_t getNumGapInsertions() const;
   hal_size_t getNumGapDeletions() const;
   hal_size_t getNumGapInsertedBases() const;
   hal_size_t getNumGapDeletedBases() const;
   bool isLast() const;
   bool isFirst() const;
   hal_index_t getLeftArrayIndex() const;
   hal_index_t getRightArrayIndex() const;
   
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

   // GappedBottomSegmentIterator Interface Methods
   GappedBottomSegmentIteratorPtr copy();
   GappedBottomSegmentIteratorConstPtr copy() const;
   void copy(GappedBottomSegmentIteratorConstPtr bs) const;
   void toParent(GappedTopSegmentIteratorConstPtr ts) const;
   bool equals(GappedBottomSegmentIteratorConstPtr other) const;
   bool hasChild() const;
   bool getChildReversed() const;
   BottomSegmentIteratorConstPtr getLeft() const;
   BottomSegmentIteratorConstPtr getRight() const;

private:
   
   bool compatible(BottomSegmentIteratorConstPtr left,
                   BottomSegmentIteratorConstPtr right) const;

   void extendRight() const;
   void extendLeft() const;

   void toLeftNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toRightNextUngapped(BottomSegmentIteratorConstPtr bs) const;
   void toLeftNextUngapped(TopSegmentIteratorConstPtr ts) const;
   void toRightNextUngapped(TopSegmentIteratorConstPtr ts) const;
   
   // keep convention of other iterators where const-ness only applies
   // to the database and not the iterator...
   mutable BottomSegmentIteratorConstPtr _left;
   mutable BottomSegmentIteratorConstPtr _right;
   mutable TopSegmentIteratorConstPtr _leftChild;
   mutable TopSegmentIteratorConstPtr _rightChild;
   mutable BottomSegmentIteratorConstPtr _temp;
   mutable BottomSegmentIteratorConstPtr _temp2;
   mutable hal_size_t _childIndex;
   mutable hal_size_t _gapThreshold;            
  
};


}
#endif
