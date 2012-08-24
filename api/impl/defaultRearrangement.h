/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDEFAULTREARRANGEMENT_H
#define _HALDEFAULTREARRANGEMENT_H

#include "halRearrangement.h"
#include "halGappedTopSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"

namespace hal {

/** 
 * Default (implementation independent) implementation of Rearrangement
 * interface 
 */
class DefaultRearrangement : public Rearrangement
{
public:
   static const hal_size_t DefaultGapThreshold;

   DefaultRearrangement(const Genome* childGenome, 
                        hal_size_t gapThreshold = DefaultGapThreshold);

   ~DefaultRearrangement();
   
   // Rearrangement Interface Methods
   ID getID() const;
   hal_size_t getLength() const;
   hal_size_t getNumContainedGaps() const;
   hal_size_t getNumContainedGapBases() const;
   TopSegmentIteratorConstPtr getLeftBreakpoint() const;
   TopSegmentIteratorConstPtr getRightBreakpoint() const;
   bool identifyFromLeftBreakpoint(TopSegmentIteratorConstPtr topSegment);
   bool identifyNext();
   hal_size_t getGapLengthThreshold() const;
   void setGapLengthThreshold(hal_size_t threshold);

private:
   
   void resetStatus(TopSegmentIteratorConstPtr topSegment);

   bool scanInversionCycle(TopSegmentIteratorConstPtr topSegment);
   bool scanInsertionCycle(TopSegmentIteratorConstPtr topSegment);
   bool scanDeletionCycle(TopSegmentIteratorConstPtr topSegment);
   bool scanTranslocationCycle(TopSegmentIteratorConstPtr topSegment);
   bool scanDuplicationCycle(TopSegmentIteratorConstPtr topSegment);

private:
   
   hal_size_t _gapThreshold;

   GappedTopSegmentIteratorConstPtr _left, _cur, _right;
   GappedBottomSegmentIteratorConstPtr _leftParent, _rightParent, _tempParent;
   TopSegmentIteratorConstPtr _top;
   
   hal_size_t _childIndex;
   const Genome* _genome;
   const Genome* _parent;

   ID _id;  
};


}
#endif
