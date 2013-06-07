/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTITERATOR_H
#define _HALSEGMENTITERATOR_H

#include "halDefs.h"
#include "halSlicedSegment.h"

namespace hal {

/** 
 * Interface for general segment iterator.  Common functionality
 * for top and bottom iterators.  A segment iterator implements 
 * the segment as well as the iterator interface. 
 */
class SegmentIterator : public virtual SlicedSegment
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

   /** move iterator to position of dna site in *genome*
    * @param position index of site in genome
    * @param slice if true, the returned iterator is sliced to correspond
    * to just the single base.  if false, the iterator corresponds to the
    * entire segment. 
    * NOTE*** this function requires up to log2(N) time in current hdf5 imp.
    * though it should be faster on average*/
   virtual void toSite(hal_index_t position, bool slice = true) const = 0;


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
