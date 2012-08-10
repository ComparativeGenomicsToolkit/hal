/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALDEFAULTREARRANGEMENT_H
#define _HALDEFAULTREARRANGEMENT_H

#include "halRearrangement.h"

namespace hal {

/** 
 * Default (implementation independent) implementation of Rearrangement
 * interface 
 */
class DefaultRearrangement : public Rearrangement
{
public:

   DefaultRearrangement();

   ~DefaultRearrangement();
   
   // Rearrangement Interface Methods
   ID getID() const;
   hal_size_t getLength() const;
   hal_size_t getNumContainedGaps() const;
   hal_size_t getNumContainedGapBases() const;
   hal_size_t getDuplicationDegree() const;
   TopSegmentIteratorConstPtr getLeftBreakpoint() const;
   TopSegmentIteratorConstPtr getRightBreakpoint() const;
   bool identifyFromLeftBreakpoint(TopSegmentIteratorConstPtr topSegment);
   bool identifyNext();
   hal_size_t getGapLengthThreshold() const;
   void setGapLengthThreshold(hal_size_t threshold);

private:

   static const hal_size_t DefaultGapThreshold;
   
   hal_size_t _gapThreshold;

   TopSegmentIteratorConstPtr _left;
   TopSegmentIteratorConstPtr _right;
   BottomSegmentIteratorConstPtr _leftParent;
   BottomSegmentIteratorConstPtr _rightParent;

   ID _id;
   hal_size_t _length;
   hal_size_t _numGaps;
   hal_size_t _numGapBases;
   hal_size_t _dupDegree;   
};


}
#endif
