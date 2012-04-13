/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENT_H
#define _HALBOTTOMSEGMENT_H

#include "halDefs.h"
#include "halSegment.h"

namespace hal {

/** 
 * Interface for a bottom segment of DNA
 */
class BottomSegment : public virtual Segment
{
public:

   /** Get the number of child genomes (note this is a number of slots
    * and that the current segment could actually have fewer children) */
   virtual hal_size_t getNumChildren() const = 0;
   
   /** Get the index of a child segment (OR NULL_INDEX if none)
    * @param i index of child to query */
   virtual hal_index_t getChildIndex(hal_size_t i) const = 0;

   /** Set the index of a child segment (OR NULL_INDEX if none)
    * @param i index of child to set 
    * @param childIndex index of segment in child to set */
   virtual void setChildIndex(hal_size_t i, hal_index_t childIndex) = 0;

   /** Get whether descent segment for ith child is mapped to the 
    * reverse complement of this segment 
    * @param i index of child to query */
   virtual hal_bool_t getChildReversed(hal_size_t i) const = 0;

   /** Set whether descent segment for ith child is mapped to the 
    * reverse complement of this segment 
    * @param i index of child to set 
    * @param isReverse flag */
   virtual void setChildReversed(hal_size_t child, hal_bool_t isReversed) = 0;

   /** Get index of top segment in samge genome that contains
    * this segment's start coordinate */
   virtual hal_index_t getTopParseIndex() const = 0;

   /** Set index of top segment in samge genome that contains
    * this segment's start coordinate 
    * @param parParseIndex index */
   virtual void setTopParseIndex(hal_index_t parseIndex) = 0;
   
   /** Get offset in associated top segment of start coordinate of 
    * this segment */
   virtual hal_offset_t getTopParseOffset() const = 0;

   /** Set offset in associated top segment of start coordinate of 
    * this segment 
    * @param parpArseOffset offset in parent */
   virtual void setTopParseOffset(hal_offset_t parseOffset) = 0;

protected:

   /** Destructor */
   virtual ~BottomSegment() = 0;
};

inline BottomSegment::~BottomSegment() {}
}
#endif
