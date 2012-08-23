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
    * the segment iterator is sliced to begin at this position 
    * (if the iterator is reversed, it moves in opposite direction) */
   virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX) const = 0;

   /** move iterator one segment to the right 
    * @param rightCutoff If the right segment contains the 
    * genome position (DNA coordinate) specified by rightCutoff
    * the segment iterator is sliced to end at this position 
    * (if the iterator is reversed, it moves in opposite direction) */
   virtual void toRight(hal_index_t rightCutoff = NULL_INDEX) const = 0;

   /** switch to segment's reverse complement */
   virtual void toReverse() const = 0;

   /** move iterator to position of dna site in *genome*
    * @param position index of site in genome
    * @param slice if true, the returned iterator is sliced to correspond
    * to just the single base.  if false, the iterator corresponds to the
    * entire segment. 
    * NOTE*** this function requires up to log2(N) time in current hdf5 imp.
    * though it should be faster on average*/
   virtual void toSite(hal_index_t position, bool slice = true) const = 0;

   /** check if there is a next paralogy */
   virtual bool hasNextParalogy() const = 0;

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
    * start position plus the startoffset.  Note that the start position
    * is currently always in FORWARD GENOME coordinates.  That said,
    * if the segment iterator is in reverse orientation, the segment
    * is read from right (startposition) to left (startposition - length -1)*/
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

   /** Is the first segment of a sequence 
    * ie has index 0 in sequence if not reversed or index n-1 if is
    * reversed */
   virtual bool isFirst() const = 0;

   /** Is the last segment of a sequence 
    * ie has index n-1 in sequence if not reversed or index 0 if is
    * reversed */
   virtual bool isLast() const = 0;

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
