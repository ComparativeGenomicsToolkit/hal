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
class DefaultGappedTopSegmentIterator : public virtual GappedTopSegmentIterator
{
public:

   DefaultGappedTopSegmentIterator(TopSegmentIteratorConstPtr left,
                                   hal_size_t gapThreshold,
                                   bool atomic);

   virtual ~DefaultGappedTopSegmentIterator();
   
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
   virtual bool isTop() const;
   virtual hal_size_t getMappedSegments(
     std::set<MappedSegmentConstPtr>& outSegments,
     const Genome* tgtGenome,
     const std::set<const Genome*>* genomesOnPath,
     bool doDupes,
     hal_size_t minLength,
     const Genome *coalescenceLimit,
     const Genome *mrca) const;
   virtual void print(std::ostream& os) const;

   // SEGMENT ITERATOR INTERFACE
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   virtual void toReverse() const;
   virtual void toReverseInPlace() const;
   virtual void toSite(hal_index_t position, bool slice = true) const;
   virtual hal_offset_t getStartOffset() const;
   virtual hal_offset_t getEndOffset() const;
   virtual void slice(hal_offset_t startOffset,
                      hal_offset_t endOffset) const;
   virtual bool getReversed() const;

   // GAPPED SEGMENT ITERATOR INTERFACE
   virtual hal_size_t getGapThreshold() const;
   virtual bool getAtomic() const;
   virtual hal_size_t getChildIndex() const;
   virtual hal_size_t getNumSegments() const;
   virtual hal_size_t getNumGaps() const;
   virtual hal_size_t getNumGapBases() const;
   virtual hal_index_t getLeftArrayIndex() const;
   virtual hal_index_t getRightArrayIndex() const;

   // GAPPED TOP SEGMENT ITERATOR INTERFACE
   virtual GappedTopSegmentIteratorPtr copy();
   virtual GappedTopSegmentIteratorConstPtr copy() const;
   virtual void copy(GappedTopSegmentIteratorConstPtr ts) const;
   virtual void toChild(GappedBottomSegmentIteratorConstPtr bs) const;
   virtual bool equals(GappedTopSegmentIteratorConstPtr other) const;
   virtual bool adjacentTo(GappedTopSegmentIteratorConstPtr other) const;
   virtual bool hasParent() const;
   virtual bool hasNextParalogy() const;
   virtual void toNextParalogy() const;
   virtual bool getParentReversed() const;
   virtual TopSegmentIteratorConstPtr getLeft() const;
   virtual TopSegmentIteratorConstPtr getRight() const;
   virtual void setLeft(TopSegmentIteratorConstPtr ts) const;
   virtual bool isCanonicalParalog() const;

protected:
   
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
   mutable bool _atomic;
  
};


}
#endif
