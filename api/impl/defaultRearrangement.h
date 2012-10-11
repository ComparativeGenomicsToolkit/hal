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
                        hal_size_t gapThreshold = DefaultGapThreshold,
                        bool atomic = false);

   virtual ~DefaultRearrangement();
   
   // Rearrangement Interface Methods
   virtual ID getID() const;
   virtual hal_size_t getLength() const;
   virtual hal_size_t getNumContainedGaps() const;
   virtual hal_size_t getNumContainedGapBases() const;
   virtual TopSegmentIteratorConstPtr getLeftBreakpoint() const;
   virtual TopSegmentIteratorConstPtr getRightBreakpoint() const;
   virtual bool identifyFromLeftBreakpoint(TopSegmentIteratorConstPtr topSegment);
   virtual bool identifyDeletionFromLeftBreakpoint(
     TopSegmentIteratorConstPtr topSegment);
   virtual std::pair<hal_index_t, hal_index_t> getDeletedRange() const;
   virtual bool identifyInsertionFromLeftBreakpoint(
     TopSegmentIteratorConstPtr topSegment);
   virtual std::pair<hal_index_t, hal_index_t> getInsertedRange() const;
   virtual bool identifyNext();
   virtual hal_size_t getGapLengthThreshold() const;
   virtual void setGapLengthThreshold(hal_size_t threshold);
   virtual bool getAtomic() const;
   virtual void setAtomic(bool atomic);

protected:
   
   virtual void resetStatus(TopSegmentIteratorConstPtr topSegment);
   
   virtual bool scanNothingCycle(TopSegmentIteratorConstPtr topSegment);
   virtual bool scanInversionCycle(TopSegmentIteratorConstPtr topSegment);
   virtual bool scanInsertionCycle(TopSegmentIteratorConstPtr topSegment);
   virtual bool scanDeletionCycle(TopSegmentIteratorConstPtr topSegment);
   virtual bool scanTranslocationCycle(TopSegmentIteratorConstPtr topSegment);
   virtual bool scanDuplicationCycle(TopSegmentIteratorConstPtr topSegment);

protected:
   
   hal_size_t _gapThreshold;
   bool _atomic;

   GappedTopSegmentIteratorConstPtr _left, _cur, _next, _right;
   GappedBottomSegmentIteratorConstPtr _leftParent, _curParent, _nextParent,
   _rightParent;
   TopSegmentIteratorConstPtr _top;
   
   hal_size_t _childIndex;
   const Genome* _genome;
   const Genome* _parent;

   ID _id;  
};


}
#endif
