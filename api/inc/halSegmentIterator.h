/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for general segment iterator.  Common functionality
 * for top and bottom iterators.  A segment iterator implements 
 * the segment as well as the iterator interface. 
 */
class SegmentIterator 
{
public:

   /** move iterator one segment to the left 
    * @param leftCutoff If the left segment contains the 
    * genome position (DNA coordinate) specified by leftCutoff
    * the segment iterator is sliced to begin at this position */
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const = 0;

   /** move iterator one segment to the right 
    * @param rightCutoff If the right segment contains the 
    * genome position (DNA coordinate) specified by rightCutoff
    * the segment iterator is sliced to end at this position */
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const = 0;

   /** switch to segment's reverse complement */
   virtual void toReverse() const = 0;

   /** move to the next paralgous segment */
   virtual void toNextParalogy() const = 0;
   
   /** Get the iterator's start offset.  This is used when moving
    * vertically and following the parse index.  Any part of the
    * segment before the start offset is ignored by the iterator */  
   virtual hal_offset_t getStartOffset() const = 0;

   /** Get the iterator's end offset.  This is used when moving
    * vertically and following the parse index.  Any part of the
    * segment after the end offset is ignored by the iterator */  
   virtual hal_offset_t getEndOffset() const = 0;

   /** Set the iterator's start and end offsets
    * @param startOffset offset from beginning of segment
    * @param endOffset offset from end of segment */
   virtual void slice(hal_offset_t startOffset = 0,
                      hal_offset_t endOffset = 0) const = 0;

   /** Get the start position of the iterator's segment.  The parse
    * offset is taken into account, so it will return the segment's
    * start position plus the startoffset */
   virtual hal_index_t getStartPosition() const = 0;

   /** Get the length of the iterator's segment.  The parse offset's
    * are taken into account.  So this will return the segment's 
    * length minues the sum of the two offsets */
   virtual hal_size_t getLength() const = 0;

   /** Check whether iterator is on segment's reverse complement */
   virtual hal_bool_t getReversed() const = 0;
   
   /** Get the DNA string corresponding to the iterator from the genome 
    * @param outString string into which the results are copied */
   virtual void getString(std::string& outString) const = 0;

   /** Determine if current segment is to the left of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool leftOf(hal_index_t genomePos) const = 0;

   /** Determine if current segment is to the right of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool rightOf(hal_index_t genomePos) const = 0;

   /** Determine if current segment is to the right of a genome coordinate
    * @param genomePos Index of DNA character in genome */
   virtual bool overlaps(hal_index_t genomePos) const = 0;

protected:
   friend class counted_ptr<SegmentIterator>;
   friend class counted_ptr<const SegmentIterator>;
   virtual ~SegmentIterator() = 0;
};

inline SegmentIterator::~SegmentIterator() {}

inline bool operator<(SegmentIteratorConstPtr segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt->leftOf(genomePos);
}

inline bool operator>(SegmentIteratorConstPtr segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt->rightOf(genomePos);
}

}
#endif
