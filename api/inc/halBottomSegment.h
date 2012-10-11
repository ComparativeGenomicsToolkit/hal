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
class BottomSegment : virtual public Segment
{
public:

   /** Get the number of child genomes (note this is a number of slots
    * and that the current segment could actually have fewer children) */
   virtual hal_size_t getNumChildren() const = 0;
   
   /** Get the index of a child segment (OR NULL_INDEX if none)
    * @param i index of child to query */
   virtual hal_index_t getChildIndex(hal_size_t i) const = 0;

   /** Get the index of a child segment (OR NULL_INDEX if none)
    * @param childGenome genome of child to query */
   virtual hal_index_t getChildIndexG(const Genome* childGenome) const = 0;

   /** Test if child segment exists 
    * @param child index of child genome */
   virtual bool hasChild(hal_size_t child) const = 0;

   /** Test if child segment exists 
    * @param child Child genome */
   virtual bool hasChildG(const Genome* childGenome) const = 0;

   /** Set the index of a child segment (OR NULL_INDEX if none)
    * @param i index of child to set 
    * @param childIndex index of segment in child to set */
   virtual void setChildIndex(hal_size_t i, hal_index_t childIndex) = 0;

   /** Get whether descent segment for ith child is mapped to the 
    * reverse complement of this segment 
    * @param i index of child to query */
   virtual bool getChildReversed(hal_size_t i) const = 0;

   /** Set whether descent segment for ith child is mapped to the 
    * reverse complement of this segment 
    * @param i index of child to set 
    * @param isReverse flag */
   virtual void setChildReversed(hal_size_t child, bool isReversed) = 0;

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

   /** Test if there is a top parse segment in the genome that contains
    * the start positon of the current iterator (should be true if not
    * in root */
   virtual bool hasParseUp() const = 0;

   /** Get the index of the child of the left neighbour of this segment
    * in the genome (use isLeft first to check if the left neighbour
    * is in the same sequence)
    * @param i index of child genome to query */
   virtual hal_index_t getLeftChildIndex(hal_size_t i) const = 0;

   /** Get the right of the child of the left neighbour of this segment
    * in the genome (use isRight first to check if the right neighbour
    * is in the same sequence)
    * @param i index of child genome to query */
   virtual hal_index_t getRightChildIndex(hal_size_t i) const = 0;

protected:
   friend class counted_ptr<BottomSegment>;
   friend class counted_ptr<const BottomSegment>;
   virtual ~BottomSegment() = 0;
};

inline BottomSegment::~BottomSegment() {}
}
#endif
