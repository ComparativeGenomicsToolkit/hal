/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALREARRANGEMENT_H
#define _HALREARRANGEMENT_H

#include <vector>
#include <string>
#include "halDefs.h"
#include "halTopSegmentIterator.h"

namespace hal {

/** 
 * Interface for a rearrangement operation.  Rearrangements are specifed by
 * pairs of breakpoints in a top segment array.  For now they are two-break
 * (balance or not) rearrangements and translocations.  Anything else will
 * will be labled complex for now. 
 */
class Rearrangement 
{
public:

   enum ID { Insertion = 0, Deletion, Duplication, Translocation, 
             Transposition, Inversion, Complex, Invalid, Gap, Nothing};

   /** Get the ID of the rearrangement */
   virtual ID getID() const = 0;
   
   /** Get the length of a rearrangement (ie number of bases that were
    * inserted, deleted, moved, etc */
   virtual hal_size_t getLength() const = 0;
   
   /** Get the number of gaps within the arrangement.  Gaps are dependent
    * on the gap threshold that can be specified below... */
   virtual hal_size_t getNumContainedGaps() const = 0;

   /** Get the total size (in bases) of gaps within the arrangement.  
    * Gaps are dependent
    * on the gap threshold that can be specified below...*/
   virtual hal_size_t getNumContainedGapBases() const = 0;

   /** Left breakpoint is specified by the right cap of the returned segment */
   virtual TopSegmentIteratorConstPtr getLeftBreakpoint() const = 0;

   /** Right breakpoint is specified by the left cap of the returned segment */
   virtual TopSegmentIteratorConstPtr getRightBreakpoint() const = 0;

   /** Identify the rearrangement by scanning along the genome, skipping 
    * gaps until the matching breeakpoint is found.  If no matching breakpint
    * is found, the rearrangement cannot be identified and false is returned */
   virtual bool identifyFromLeftBreakpoint(TopSegmentIteratorConstPtr
                                           topSegment) = 0;
   
   /** Start scanning from last identified breakpoint */
   virtual bool identifyNext() = 0;

   /** Return the maximum size of a simple indel such that it will be 
    * considered a gap (and factored out of breakpoint computations) */
   virtual hal_size_t getGapLengthThreshold() const = 0;

   /** Specify the maximum size of a simple indel such that it will be 
    * considered a gap (and factored out of breakpoint computations) */
   virtual void setGapLengthThreshold(hal_size_t threshold) = 0;

protected:
   friend class counted_ptr<Rearrangement>;
   friend class counted_ptr<const Rearrangement>;
   /** Destructor */
   virtual ~Rearrangement() = 0;
};

inline Rearrangement::~Rearrangement() {}

}
#endif
