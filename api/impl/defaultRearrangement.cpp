/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <deque>
#include "defaultRearrangement.h"
#include "hal.h"

using namespace std;
using namespace hal;

// maximum size a simple indel can be to be considered a gap (and not
// a rearrangement)
const hal_size_t DefaultRearrangement::DefaultGapThreshold = 10;

DefaultRearrangement::DefaultRearrangement(const Genome* childGenome,
                                           hal_size_t gapThreshold) :
  _genome(childGenome),
  _parent(NULL),
  _gapThreshold(gapThreshold),
  _childIndex(1000),
  _id(Invalid),
  _length(0),
  _numGaps(0),
  _numGapBases(0),
  _dupDegree(0)
{
  _parent = childGenome->getParent();
  assert(_parent != NULL);
  // just allocating here -- need to be properly init
  _left = _genome->getGappedTopSegmentIterator(0, _gapThreshold);
  _right = _left->copy();
  _leftParent = _parent->getGappedBottomSegmentIterator(0, 0, _gapThreshold);
  _rightParent = _leftParent->copy();
}

DefaultRearrangement::~DefaultRearrangement()
{

}
   
// Rearrangement Interface Methods
DefaultRearrangement::
ID DefaultRearrangement::getID() const
{
  return _id;
}

hal_size_t DefaultRearrangement::getLength() const
{
  return _length;
}

hal_size_t DefaultRearrangement::getNumContainedGaps() const
{
  return _numGaps;
}

hal_size_t DefaultRearrangement::getNumContainedGapBases() const
{
  return _numGapBases;
}

hal_size_t DefaultRearrangement::getDuplicationDegree() const
{
  return _dupDegree;
}

TopSegmentIteratorConstPtr DefaultRearrangement::getLeftBreakpoint() const
{
  assert(_left->getReversed() == false);
  return _left->getLeft();
}

TopSegmentIteratorConstPtr DefaultRearrangement::getRightBreakpoint() const
{
  assert(_left->getReversed() == false);
  return _right->getRight();
}

bool DefaultRearrangement::identifyFromLeftBreakpoint(
  TopSegmentIteratorConstPtr topSegment)
{
  if (scanInversionCycle(topSegment) == true)
  {
    _id = Inversion;
  }
  else if (scanInsertionCycle(topSegment) == true)
  {
    _id = Insertion;
  }
  else if (scanDeletionCycle(topSegment) == true)
  {
    _id = Deletion;
  }
  else
  {
    _id = Complex;
  }
  
  return true;
}

bool DefaultRearrangement::identifyNext()
{
  if (_right->isLast())
  {
    return false;
  }
  _right->toRight();
  _left->copy(_right);
  assert(_left->getReversed() == false);
  bool res = identifyFromLeftBreakpoint(_left->getLeft());
  return res;
}

hal_size_t DefaultRearrangement::getGapLengthThreshold() const
{
  return _gapThreshold;
}

void DefaultRearrangement::setGapLengthThreshold(hal_size_t threshold)
{
  _gapThreshold = threshold;
}

void DefaultRearrangement::resetStatus(TopSegmentIteratorConstPtr topSegment)
{  
  _id = Invalid;
  _length = 0;
  _numGaps = 0;
  _numGapBases = 0;
  _dupDegree = 0;
  assert(topSegment.get());
  _genome = topSegment->getTopSegment()->getGenome();
  _parent = _genome->getParent();
  assert(_parent != NULL);

  _left->setLeft(topSegment);
  _right->copy(_left);
}




// Segment is an inverted descendant of another Segment
// (Note that it can still be part of a more complex operation 
// such as a transposition -- which would have to be tested for
// in the insertion cycle code.  so this method is insufficient to
// return an ID=Inversion)
bool DefaultRearrangement::scanInversionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  return _left->getParentReversed();
}

bool DefaultRearrangement::scanInsertionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  return false;
}

bool DefaultRearrangement::scanDeletionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  return false;
}
