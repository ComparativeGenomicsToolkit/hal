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
   hal_size_t getDuplicationDegree() const;
   TopSegmentIteratorConstPtr getLeftBreakpoint() const;
   TopSegmentIteratorConstPtr getRightBreakpoint() const;
   bool identifyFromLeftBreakpoint(TopSegmentIteratorConstPtr topSegment);
   bool identifyNext();
   hal_size_t getGapLengthThreshold() const;
   void setGapLengthThreshold(hal_size_t threshold);

private:
   
   void resetStatus(TopSegmentIteratorConstPtr topSegment);

   bool isGap(TopSegmentIteratorConstPtr topSegment);
   bool isGap(BottomSegmentIteratorConstPtr bottomSegment);

   bool isForwardAligned(TopSegmentIteratorConstPtr inLeft, 
                         TopSegmentIteratorConstPtr inRight,
                         BottomSegmentIteratorConstPtr inParentLeft,
                         BottomSegmentIteratorConstPtr inParentRight);

   bool isReverseAligned(TopSegmentIteratorConstPtr inLeft, 
                        TopSegmentIteratorConstPtr inRight,
                        BottomSegmentIteratorConstPtr inParentLeft,
                        BottomSegmentIteratorConstPtr inParentRight);

   // traversal shortcuts used to skip gaps.
   // left is always left in genome coordiantes (even if iterator is 
   // inverted)
   void findRight(TopSegmentIteratorConstPtr inLeft,
                  TopSegmentIteratorConstPtr& outRight);
   void findLeft(TopSegmentIteratorConstPtr inRight,
                 TopSegmentIteratorConstPtr& outLeft);
   void findRight(BottomSegmentIteratorConstPtr inLeft,
                  BottomSegmentIteratorConstPtr& outRight);
   void findLeft(BottomSegmentIteratorConstPtr inRight,
                 BottomSegmentIteratorConstPtr& outLeft);

   void findParents(TopSegmentIteratorConstPtr inLeft, 
                    TopSegmentIteratorConstPtr inRight,
                    BottomSegmentIteratorConstPtr& outParentLeft,
                    BottomSegmentIteratorConstPtr& outParentRight);

   void findChilds(BottomSegmentIteratorConstPtr inParentLeft, 
                   BottomSegmentIteratorConstPtr inParentRight,
                   TopSegmentIteratorConstPtr& outLeft,
                   TopSegmentIteratorConstPtr& outRight);

   // Inversion Module 
   bool scanInversionCycle(TopSegmentIteratorConstPtr topSegment);

   // Insertion Modules
   // see defaultRearrangementIns.cpp
   bool scanInsertionCycle(TopSegmentIteratorConstPtr topSegment);
   bool middleInsertion();
   bool startInsertion();
   bool endInsertion();
   
   // Deletion Modules
   // see defaultRearrangementDel.cpp
   bool scanDeletionCycle(TopSegmentIteratorConstPtr topSegment);
   bool middleDeletion();
   bool startDeletion();
   bool endDeletion();

private:
   
   hal_size_t _gapThreshold;

   GappedTopSegmentIteratorConstPtr _left, _right;
   GappedBottomSegmentIteratorConstPtr _leftParent, _rightParent;

   TopSegmentIteratorConstPtr _startLeft, _startRight;
   TopSegmentIteratorConstPtr _startLeft2, _startRight2;
   TopSegmentIteratorConstPtr _endLeft, _endRight;
   BottomSegmentIteratorConstPtr _parentStartLeft, _parentStartRight;
   BottomSegmentIteratorConstPtr _parentEndLeft, _parentEndRight;
   
   hal_size_t _childIndex;
   const Genome* _genome;
   const Genome* _parent;

   ID _id;
   hal_size_t _length;
   hal_size_t _numGaps;
   hal_size_t _numGapBases;
   hal_size_t _dupDegree;   
  
};


}
#endif
