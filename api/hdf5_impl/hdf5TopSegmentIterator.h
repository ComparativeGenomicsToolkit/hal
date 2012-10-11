/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5TOPSEGMENTITERATOR_H
#define _HDF5TOPSEGMENTITERATOR_H

#include <H5Cpp.h>
#include "halTopSegmentIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5TopSegment.h"

namespace hal {

class HDF5BottomSegmentIterator;

class HDF5TopSegmentIterator : public TopSegmentIterator
{
   friend class HDF5BottomSegmentIterator;
public:
   
   HDF5TopSegmentIterator(HDF5Genome* genome, hal_index_t index,
                          hal_offset_t startOffset = 0, 
                          hal_offset_t endOffset = 0,
                          bool inverted = false);
   ~HDF5TopSegmentIterator();
   
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

    // TOP SEGMENT INTERFACE
   hal_index_t getParentIndex() const;
   bool hasParent() const;
   void setParentIndex(hal_index_t parIdx);
   bool getParentReversed() const;
   void setParentReversed(bool isReversed);
   hal_index_t getBottomParseIndex() const;
   void setBottomParseIndex(hal_index_t botParseIdx);
   hal_offset_t getBottomParseOffset() const;
   bool hasParseDown() const;
   hal_index_t getNextParalogyIndex() const;
   bool hasNextParalogy() const;
   void setNextParalogyIndex(hal_index_t parIdx);
   hal_index_t getLeftParentIndex() const;
   hal_index_t getRightParentIndex() const;

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

   // TOP SEGMENT ITERATOR INTERFACE 
   TopSegmentIteratorPtr copy();
   TopSegmentIteratorConstPtr copy() const;
   void copy(TopSegmentIteratorConstPtr ts) const;
   void toChild(BottomSegmentIteratorConstPtr bs, 
                        hal_size_t child) const;
   void toChildG(BottomSegmentIteratorConstPtr bs, 
                 const Genome* childGenome) const;
   void toParseUp(BottomSegmentIteratorConstPtr bs) const;
   TopSegment* getTopSegment();
   const TopSegment* getTopSegment() const;
   bool equals(TopSegmentIteratorConstPtr other) const;
   void toNextParalogy() const;

private:
   bool inRange() const;
   
private:
   HDF5TopSegment _topSegment;
   mutable hal_offset_t _startOffset;
   mutable hal_offset_t _endOffset;
   mutable bool _reversed;
};

inline bool HDF5TopSegmentIterator::inRange() const
{
  return _topSegment._index >= 0 && _topSegment._index < 
     (hal_index_t)_topSegment._genome->getNumTopSegments();
}

}
#endif
