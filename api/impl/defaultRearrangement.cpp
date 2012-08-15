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

DefaultRearrangement::DefaultRearrangement() :
  _gapThreshold(DefaultGapThreshold),
  _childIndex(1000),
  _genome(NULL),
  _parent(NULL),
  _id(Invalid),
  _length(0),
  _numGaps(0),
  _numGapBases(0),
  _dupDegree(0)
{

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
  return _startLeft;
}

TopSegmentIteratorConstPtr DefaultRearrangement::getRightBreakpoint() const
{
  return _endLeft;
}

bool DefaultRearrangement::identifyFromLeftBreakpoint(
  TopSegmentIteratorConstPtr topSegment)
{
  // initalise the object
  resetStatus(topSegment);
  
  // Determine cycle by 1) scanning right (skipping gaps).
  // Scan up (first segment with a parent)
  // Scan left (skipping gaps)
  // scan down (first segment with a child)
  
  return false;
}

bool DefaultRearrangement::identifyNext()
{
  _startLeft->copy(_endRight);
  _startLeft->toRight();
  return identifyFromLeftBreakpoint(_startLeft);
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
  assert(topSegment.get() != NULL);
  if (_startLeft.get() == NULL)
  {
    _startLeft = topSegment->copy();
  }
  else
  {
    _startLeft->copy(topSegment);
  }
  if (_startLeft->getReversed() == true)
  {
    _startLeft->toReverse();
  }

  // only reinitialize all the genome-related machinery if the 
  // input genome is different than the current genome
  if (_genome != _startLeft->getTopSegment()->getGenome())
  {
    _genome = _startLeft->getTopSegment()->getGenome();
    _parent = _genome->getParent();
    assert(_parent != NULL);
    hal_size_t numChildren = _parent->getNumChildren();
    _childIndex = numChildren;
    for (hal_size_t i = 0; i < numChildren; ++i)
    {
      if (_parent->getChild(i) == _genome)
      {
        _childIndex = i;
      }
    }
    assert(_childIndex < numChildren);
    _startLeft = _genome->getTopSegmentIterator();
    _startLeft2 = _genome->getTopSegmentIterator();
    _startRight = _genome->getTopSegmentIterator();
    _startRight2 = _genome->getTopSegmentIterator();
    _parentStartLeft = _parent->getBottomSegmentIterator();
    _parentEndLeft = _parent->getBottomSegmentIterator();
    _parentStartRight = _parent->getBottomSegmentIterator();
    _parentEndRight = _parent->getBottomSegmentIterator();
  }
  findRight(_startLeft, _startRight);

  assert(_genome != NULL);
  assert(_startLeft.get() != NULL);
  assert(_parent != NULL);
  assert(_parent->getChild(_childIndex) == _genome);
}

void DefaultRearrangement::findRight(TopSegmentIteratorConstPtr inLeft,
                                     TopSegmentIteratorConstPtr outRight)
{
  outRight->copy(inLeft);
  if (outRight->getReversed() == true)
  {
    outRight->toReverse();
  }

  while (outRight->getTopSegment()->isLast() == false)
  {
    outRight->toRight();
    if (outRight->getTopSegment()->isGapInsertion() == false &&
        outRight->getTopSegment()->isSimpleInversion() == false)
    {
      break;
    }
  }
  assert(outRight->getTopSegment()->getArrayIndex() >= 
         inLeft->getTopSegment()->getArrayIndex());
}

void DefaultRearrangement::findLeft(TopSegmentIteratorConstPtr inRight,
                                    TopSegmentIteratorConstPtr outLeft)
{
  outLeft->copy(inRight);
  if (outLeft->getReversed() == true)
  {
    outLeft->toReverse();
  }

  while (outLeft->getTopSegment()->isFirst() == false)
  {
    outLeft->toLeft();
    if (outLeft->getTopSegment()->isGapInsertion() == false &&
        outLeft->getTopSegment()->isSimpleInversion() == false)
    {
      break;
    }
  }
  assert(inRight->getTopSegment()->getArrayIndex() >= 
         outLeft->getTopSegment()->getArrayIndex());
}

void DefaultRearrangement::findRight(BottomSegmentIteratorConstPtr inLeft,
                                     BottomSegmentIteratorConstPtr outRight)
{
  outRight->copy(inLeft);
  if (outRight->getReversed() == true)
  {
    outRight->toReverse();
  }

  while (outRight->getBottomSegment()->isLast() == false)
  {
    outRight->toRight();
    if (outRight->getBottomSegment()->isGapDeletion(_childIndex) == false &&
        outRight->getBottomSegment()->isSimpleInversion(_childIndex) == false)
    {
      outRight->toLeft();
      break;
    }
    else
    {
      outRight->toRight();
    }
  }
  assert(outRight->getBottomSegment()->getArrayIndex() >= 
         inLeft->getBottomSegment()->getArrayIndex());
}

void DefaultRearrangement::findLeft(BottomSegmentIteratorConstPtr inRight,
                                    BottomSegmentIteratorConstPtr outLeft)
{
  outLeft->copy(inRight);
  if (outLeft->getReversed() == true)
  {
    outLeft->toReverse();
  }

  while (outLeft->getBottomSegment()->isFirst() == false)
  {
    outLeft->toLeft();
    if (outLeft->getBottomSegment()->isGapDeletion(_childIndex) == false &&
        outLeft->getBottomSegment()->isSimpleInversion(_childIndex) == false)
    {
      outLeft->toRight();
      break;
    }
    else
    {
      outLeft->toLeft();
    }
  }
  assert(inRight->getBottomSegment()->getArrayIndex() >= 
         outLeft->getBottomSegment()->getArrayIndex());
}

void 
DefaultRearrangement::findParents(TopSegmentIteratorConstPtr inLeft, 
                                  TopSegmentIteratorConstPtr inRight,
                                  BottomSegmentIteratorConstPtr outParentLeft,
                                  BottomSegmentIteratorConstPtr outParentRight)
{
  assert(inLeft->hasParent() && inRight->hasParent());
  assert(inLeft->getTopSegment()->getParentReversed() == 
         inRight->getTopSegment()->getParentReversed());

  outParentLeft->toParent(inLeft);
  outParentRight->toParent(inRight);
  if (outParentLeft->getBottomSegment()->getArrayIndex() > 
      outParentRight->getBottomSegment()->getArrayIndex())
  {
    assert(inLeft->getTopSegment()->getParentReversed() == true);
    swap(outParentLeft, outParentRight);
  }
  if (outParentLeft->getReversed())
  {
    outParentLeft->toReverse();
    outParentRight->toReverse();
  }
  assert(outParentLeft->getReversed() == outParentRight->getReversed());


#ifndef NDEBUG
  TopSegmentIteratorConstPtr t = inLeft->copy();
  BottomSegmentIteratorConstPtr b = outParentLeft->copy();
  findRight(inLeft, t);
  assert(t->equals(inRight));
  findLeft(inRight, t);
  assert(t->equals(inLeft));
  findRight(outParentLeft, b);
  assert(b->equals(outParentRight));
  findLeft(outParentRight, b);
  assert(b->equals(outParentLeft));
#endif
  
}

void 
DefaultRearrangement::findChilds(BottomSegmentIteratorConstPtr inParentLeft, 
                                 BottomSegmentIteratorConstPtr inParentRight,
                                 TopSegmentIteratorConstPtr outLeft,
                                 TopSegmentIteratorConstPtr outRight)
{
  assert(inParentLeft->hasChild(_childIndex) && 
         inParentRight->hasChild(_childIndex));
  assert(inParentLeft->getBottomSegment()->getChildReversed(_childIndex) == 
         inParentRight->getBottomSegment()->getChildReversed(_childIndex));

  outLeft->toChild(inParentLeft, _childIndex);
  outRight->toChild(inParentRight, _childIndex);
  if (outLeft->getTopSegment()->getArrayIndex() > 
      outRight->getTopSegment()->getArrayIndex())
  {
    assert(inParentLeft->getBottomSegment()->getChildReversed(_childIndex) == 
           true);
    swap(outLeft, outRight);
  }
  if (outLeft->getReversed())
  {
    outLeft->toReverse();
    outRight->toReverse();
  }
  assert(outLeft->getReversed() == outRight->getReversed());

#ifndef NDEBUG
  TopSegmentIteratorConstPtr t = outLeft->copy();
  BottomSegmentIteratorConstPtr b = inParentLeft->copy();
  findRight(outLeft, t);
  assert(t->equals(outRight));
  findLeft(outRight, t);
  assert(t->equals(outLeft));
  findRight(inParentLeft, b);
  assert(b->equals(inParentRight));
  findLeft(inParentRight, b);
  assert(b->equals(inParentLeft));
#endif

}

// try to find a cycle like the following (barring gaps)
// right one segment
// up one segment
// left one segment
// down one segment
bool DefaultRearrangement::scanInversionCycle()
{
  // move right one segment
  if (_startLeft->getTopSegment()->isLast() == true)
  {
    return false;
  }
  assert(_startLeft->getReversed() == false);
  findRight(_startLeft, _startRight);
  _endLeft->copy(_startRight);
  _endLeft->toRight();
  findRight(_endLeft, _endRight);

  // move up one segment
  if (_endLeft->hasParent() == false || 
      _endLeft->getTopSegment()->getParentReversed() == false)
  {
    return false;
  }
  findParents(_endLeft, _endRight, _parentEndLeft, _parentEndRight);
  
  // move left one segment
  _parentStartLeft->copy(_parentEndLeft);
  _parentStartRight->copy(_parentEndRight);
  if (_parentEndLeft->getBottomSegment()->isFirst() == true)
  {
    return true;
  }
  _parentStartLeft->toLeft();
  findLeft(_parentStartLeft, _parentStartRight);
  swap(_parentStartLeft, _parentStartRight);
  
  // move down one segment
  if (_parentStartLeft->hasChild(_childIndex) == true)
  {
    findChilds(_parentStartLeft, _parentStartRight, _startLeft2, _startRight2);
  }
  
  if (_startLeft->equals(_startLeft2))
  {
    assert(_startRight->equals(_startRight2));
    return true;
  }

  return false;
}

// try to find a cycle like the following (barring gaps)
// up one segment
// right one segment
// down one segment
// left two segments
bool DefaultRearrangement::scanInsertionCycle()
{
  return false;
}
