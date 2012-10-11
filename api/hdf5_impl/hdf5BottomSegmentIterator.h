/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5BOTTOMSEGMENTITERATOR_H
#define _HDF5BOTTOMSEGMENTITERATOR_H

#include <limits>
#include <H5Cpp.h>
#include "halBottomSegmentIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5BottomSegment.h"

namespace hal {

class HDF5TopSegmentIterator;

class HDF5BottomSegmentIterator : public BottomSegmentIterator
{

   friend class HDF5TopSegmentIterator;

public:
   
   HDF5BottomSegmentIterator(HDF5Genome* genome, hal_index_t index,
                             hal_size_t startOffset = 0, 
                             hal_size_t endOffset = 0,
                             bool inverted = false);
   ~HDF5BottomSegmentIterator();
   
   // SEGMENT INTERFACE
   const Genome* getGenome() const;
   Genome* getGenome();
   const Sequence* getSequence() const;
   Sequence* getSequence();
   hal_index_t getStartPosition() const;
   hal_index_t getEndPosition() const;
   hal_size_t getLength() const;
   void getString(std::string& outString) const;
   void setCoordinates(hal_index_t startPos, hal_size_t length);
   hal_index_t getArrayIndex() const;
   bool leftOf(hal_index_t genomePos) const;
   bool rightOf(hal_index_t genomePos) const;
   bool overlaps(hal_index_t genomePos) const;
   bool isFirst() const;
   bool isLast() const;
   
   // BOTTOM SEGMENT INTERFACE
   hal_size_t getNumChildren() const;
   hal_index_t getChildIndex(hal_size_t i) const;
   hal_index_t getChildIndexG(const Genome* childGenome) const;
   bool hasChild(hal_size_t child) const;
   bool hasChildG(const Genome* childGenome) const;
   void setChildIndex(hal_size_t i, hal_index_t childIndex);
   bool getChildReversed(hal_size_t i) const;
   void setChildReversed(hal_size_t child, bool isReversed);
   hal_index_t getTopParseIndex() const;
   void setTopParseIndex(hal_index_t parseIndex);
   hal_offset_t getTopParseOffset() const;
   bool hasParseUp() const;
   hal_index_t getLeftChildIndex(hal_size_t i) const;
   hal_index_t getRightChildIndex(hal_size_t i) const;

   // SEGMENT ITERATOR INTERFACE 
   void toLeft(hal_index_t leftCutoff = NULL_INDEX) const;
   void toRight(hal_index_t rightCutoff = NULL_INDEX) const;
   void toReverse() const;
   void toSite(hal_index_t position, bool slice = true) const;
   hal_offset_t getStartOffset() const;
   hal_offset_t getEndOffset() const;
   void slice(hal_offset_t startOffset ,
                      hal_offset_t endOffset ) const;
   bool getReversed() const;

   // BOTTOM SEGMENT ITERATOR INTERFACE
   BottomSegmentIteratorPtr copy();
   BottomSegmentIteratorConstPtr copy() const;
   void copy(BottomSegmentIteratorConstPtr bs) const;
   void toParent(TopSegmentIteratorConstPtr ts) const; 
   void toParseDown(TopSegmentIteratorConstPtr ts) const;
   BottomSegment* getBottomSegment();
   const BottomSegment* getBottomSegment() const;
   bool equals(BottomSegmentIteratorConstPtr other) const;

private:
   bool inRange() const;

private:
   HDF5BottomSegment _bottomSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable bool _reversed;
};

inline bool HDF5BottomSegmentIterator::inRange() const
{
  return _bottomSegment._index >= 0 && _bottomSegment._index < 
     (hal_index_t)_bottomSegment._genome->getNumBottomSegments();
}

}
#endif
