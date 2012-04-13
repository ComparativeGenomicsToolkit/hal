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
class TopSegment : public Segment
{
public:
  
   /** Get index of the homologous segmenet in the ancestral genome */
   virtual hal_index_t getParentIndex() const = 0;
   
   /** Set the index of the homologous segment in the ancestra genome 
    * @param parIdx parent index to set */
   virtual void setParentIndex(hal_index_t parIdx) = 0;

   /** Check whether segment is mapped to parent's reverse complement */
   virtual hal_bool_t getParentReversed() const = 0;

   /** Set whether segment is mapped to parent's reverse complement 
    * @param isReversed Flag */
   virtual void setParentReversed(hal_bool_t isReversed) = 0;

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

   /** Set the offset in the bottom parse segment that aligns with the
    * start coordinate of this segment 
    * @param botParseOffset offset */
   virtual void setBottomParseOffset(hal_offset_t botParseOffset) = 0;

protected:

   /** Destructor */
   virtual ~TopSegment() = 0;
};

inline TopSegment::~TopSegment() {}
}
#endif
