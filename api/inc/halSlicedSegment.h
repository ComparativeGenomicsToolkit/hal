/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSLICEDSEGMENT_H
#define _HALSLICEDSEGMENT_H

#include "halDefs.h"
#include "halSegment.h"
#include "halMappedSegmentContainers.h"

namespace hal {

class MappedSegment;

/** 
 * Interface for a sliced segment.  This extends the segment interface
 * by allowing slicing (accessing just subintervals of the segment), 
 * along with reversing. 
 */
class SlicedSegment: public virtual Segment
{
public:
   /** Destructor */
    virtual ~SlicedSegment() {
    }

   /** switch to segment's reverse complement */
   virtual void toReverse() = 0;

   /** switch to segment's reverse complement without affecting the
    * coordinates in the forward strand.  Unless the segment is sliced
    * this will have an identical effect to toReverse(). Both methods
    * can be useful in different situations but the distinction is 
    * confusing.  For example, if the segment represents range [0, 10] 
    * on the forward strand but is sliced to to the subregion [3,8], 
    * then toReverse() will result in the the region [7,2], but to
    * toReverseInPlace() would yield [8,3].  */
   virtual void toReverseInPlace() = 0;

   /** Get the start offset of the slice in the segment */  
   virtual hal_offset_t getStartOffset() const = 0;

   /** Get the start offset of the slice in the segment */  
   virtual hal_offset_t getEndOffset() const = 0;

   /** Set the start and end offsets
    * @param startOffset offset from beginning of segment
    * @param endOffset offset from end of segment */
   virtual void slice(hal_offset_t startOffset = 0,
                      hal_offset_t endOffset = 0) = 0;


   /** Check whether iterator is on segment's reverse complement */
   virtual bool getReversed() const = 0;
};

inline bool operator<(const SlicedSegment& segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt.leftOf(genomePos);
}

inline bool operator>(const SlicedSegment& segmentIt,
                      hal_index_t genomePos) 
{
  return segmentIt.rightOf(genomePos);
}

}
#endif
// Local Variables:
// mode: c++
// End:
