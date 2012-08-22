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
  _cur = _genome->getGappedTopSegmentIterator(0, _gapThreshold);
  _left = _cur->copy();
  _right = _left->copy();
  _leftParent = _parent->getGappedBottomSegmentIterator(0, 0, _gapThreshold);
  _rightParent = _leftParent->copy();
  _tempParent = _leftParent->copy();
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
  assert(_cur->getReversed() == false);
  return _cur->getLeft();
}

TopSegmentIteratorConstPtr DefaultRearrangement::getRightBreakpoint() const
{
  assert(_cur->getReversed() == false);
  return _cur->getRight();
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
  if (_cur->isLast())
  {
    return false;
  }
  assert(_cur->getReversed() == false);
  // don't like this.  need to refactor interface to make better
  // use of gapped iterators
  TopSegmentIteratorConstPtr ts = _cur->getRight();
  ts->toRight();
  bool res = identifyFromLeftBreakpoint(ts);
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

  _cur->setLeft(topSegment);
  _left->copy(_cur);
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
  return _cur->getParentReversed();
}

// If true, _cur will store the insertion 'candidate'
// It must be further verified that this segment has no parent to
// distinguish between destination of transposition and insertion. 
bool DefaultRearrangement::scanInsertionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  bool first = _cur->isFirst();
  bool last = _cur->isLast();
  if (first && last)
  {
    return false;
  }

  // Case 1a) current segment is left endpoint.  we consider insertion
  // if right neighbour has parent
  if (first)
  {
    _right->toRight();
    return _right->hasParent();
  }

  // Case 1b) current segment is right endpoint.  we consider insertion
  // if left neighbour has parent
  else if (last)
  {
    _left->toLeft();
    return _left->hasParent();
  }

  // Case 2) current segment has a left neigbhour and a right neigbour
  else
  {
    _left->toLeft();
    _right->toRight();
    _leftParent->toParent(_left);
    _rightParent->toParent(_right);
    // Case 2a) Parents are adjacent
    if (_leftParent->adjacentTo(_rightParent))
    {
      return true;
    }
    // Case 2b) Left parent is endpoint
    else if (_leftParent->isFirst() || _leftParent->isLast())
    {
      return _leftParent->getSequence() == _rightParent->getSequence();
    }
    
    // Case 2c) Right parent is endpoint
    else if (_rightParent->isFirst() || _rightParent->isLast())
    {
      return _leftParent->getSequence() == _rightParent->getSequence();
    }
  }

  return false;
}

// If true, _leftParent will store the deletion 'candidate'
// It must be further verified that this segment has no child to
// distinguish between source of transposition and deletion. 
bool DefaultRearrangement::scanDeletionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  bool first = _cur->isFirst();
  bool last = _cur->isLast();
  if (_cur->hasParent() == false || (first && last))
  {
    return false;
  }

  // Case 1) current segment is a right endpoint.  we consider delection
  // if parent has neighbour
  if (last)
  {
    _leftParent->toParent(_cur);
    if (_leftParent->isFirst() == false)
    {
      _leftParent->toLeft();
      return true;
    }
    if (_leftParent->isLast() == false)
    {
      _leftParent->toRight();
      return true;
    }
  }

  // Case 2) Try to find deletion cycle by going right-up-left-left-down
  else
  {
    _leftParent->toParent(_cur);
    _right->toRight();
    if (_right->hasParent() == false)
    {
      return false;
    }
    _rightParent->toParent(_right); 
    
    if (_leftParent->getSequence() == _rightParent->getSequence())
    {
      // don't care about inversions
      // so we make sure left is left of right and they are both positive
      if (_leftParent->getReversed() == true)
      {
        _leftParent->toReverse();
      }
      if (_rightParent->getReversed() == true)
      {
        _rightParent->toReverse();
      }
      if (_rightParent->getLeftArrayIndex() < _leftParent->getLeftArrayIndex())
      {
        swap(_leftParent, _rightParent);
      }

      if (_leftParent->isLast())
      {
        return false;
      }
      
      _leftParent->toRight();
      return _leftParent->adjacentTo(_rightParent);
    }
  }

  return false;
}

// If true, _leftParent will store the swapped segment (and _cur will store)
// the other half
// NEED TO REVISE WITH STRONGER CRITERIA -- right now any operation
// next to an endpoint can get confused with a translocation.  
bool DefaultRearrangement::scanTranslocationCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  bool first = _cur->isFirst();
  bool last = _cur->isLast();
  if (_cur->hasParent() == false || (!first && !last))
  {
    return false;
  }

  _leftParent->toParent(_cur);
  bool pFirst = _leftParent->isFirst();
  bool pLast = _leftParent->isLast();
  _rightParent->copy(_leftParent);

  first ? _right->toRight() : _right->toLeft();
  pFirst ? _rightParent->toRight() : _rightParent->toLeft();

  if (_right->hasParent() == false)
  {
    return true;
  }
  else
  {
    _tempParent->toParent(_right);
    return _tempParent->equals(_rightParent);
  }
  return false;
}

// leaves duplication on _cur and _right
bool DefaultRearrangement::scanDuplicationCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  
  if (_cur->hasNextParalogy() == true)
  {
    _right->toNextParalogy();
    if (_right->getLeftArrayIndex() > _cur->getLeftArrayIndex())
    {
      return true;
    }
  }
  return false;
}
