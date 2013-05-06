/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENT_H
#define _HALTOPSEGMENT_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for a top segment of DNA
 */
class TopSegment : virtual public Segment
{
public:
  
   /** Get index of the homologous segmenet in the ancestral genome */
   virtual hal_index_t getParentIndex() const = 0;

   /* Check if segment has a parent bottom segment */
   virtual bool hasParent() const = 0;
   
   /** Set the index of the homologous segment in the ancestra genome 
    * @param parIdx parent index to set */
   virtual void setParentIndex(hal_index_t parIdx) = 0;

   /** Check whether segment is mapped to parent's reverse complement */
   virtual bool getParentReversed() const = 0;

   /** Set whether segment is mapped to parent's reverse complement 
    * @param isReversed Flag */
   virtual void setParentReversed(bool isReversed) = 0;

   /** Get the index of the bottom segment in genome that contains the
    * start coordinate of this top segment */
   virtual hal_index_t getBottomParseIndex() const = 0;

   /** Set the index of the bototm segment in the genome that contains the
    * start coordinate of this top segment 
    * @param botParseIndx index to set */
   virtual void setBottomParseIndex(hal_index_t botParseIdx) = 0;

   /** Get the offset in the bottom parse segment that aligns with the
    * start coordinate of this segment */
   virtual hal_offset_t getBottomParseOffset() const = 0;

   /* Check if segment has a "down" parse bottom segment in the same genome.
    * This segment should exist by definition unless in a leaf genome */
   virtual bool hasParseDown() const = 0;

   /** Get index of the next paralogous segment in the genome */
   virtual hal_index_t getNextParalogyIndex() const = 0;

   /* Check if top segment has a paralgy in the same genome */
   virtual bool hasNextParalogy() const = 0;
   
   /** Set index of the next paralogous segment in the genome 
    * @param parIdx of next segment in same genome that is 
    * homologous to this segment */
   virtual void setNextParalogyIndex(hal_index_t parIdx) = 0;

   /** Get the index of the parent of the left neighbour of this segment
    * in the genome (use isLeft first to check the left  
    * current segment is the first segment in a sequence) */
   virtual hal_index_t getLeftParentIndex() const = 0;

   /** Get the index of the parent of the right neighbour of this segment
    * in the genome (use isRight right to check the left  
    * current segment is the last segment in a sequence) */
   virtual hal_index_t getRightParentIndex() const = 0;

   /** For every set of paralogous top segments in a given genome, we identify a
    * unique segment as the canonical reference.  This is the segment that 
    * will be traversed when disabling duplication edges when mapping via 
    * the column iterator or mappedSegment interface.  It is also the segment
    * that is connected from its parent's down edge.*/
   virtual bool isCanonicalParalog() const = 0;

protected:
   friend class counted_ptr<TopSegment>;
   friend class counted_ptr<const TopSegment>;
   virtual ~TopSegment() = 0;
};

inline TopSegment::~TopSegment() {}
}
#endif
