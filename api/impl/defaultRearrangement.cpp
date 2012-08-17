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
  if (_endRight->getTopSegment()->getArrayIndex() == 
      _genome->getNumTopSegments() - 1)
  {
    return false;
  }
  else
  {
    _startLeft->copy(_endRight);
    _startLeft->toRight();
    bool res = identifyFromLeftBreakpoint(_startLeft);
    return res;
  }
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
    _endLeft = _genome->getTopSegmentIterator();
    _endRight = _genome->getTopSegmentIterator();
    _parentStartLeft = _parent->getBottomSegmentIterator();
    _parentEndLeft = _parent->getBottomSegmentIterator();
    _parentStartRight = _parent->getBottomSegmentIterator();
    _parentEndRight = _parent->getBottomSegmentIterator();
  }
  findRight(_startLeft, _startRight);
  _endLeft->copy(_startRight);
  _endRight->copy(_startRight);

  assert(_genome != NULL);
  assert(_startLeft.get() != NULL);
  assert(_parent != NULL);
  assert(_parent->getChild(_childIndex) == _genome);
}

bool DefaultRearrangement::isGap(TopSegmentIteratorConstPtr topSegment)
{
  return topSegment->hasParent() == false &&
     topSegment->getLength() <= _gapThreshold;
  /*
  return (topSegment->getTopSegment()->isGapInsertion() ||
     topSegment->getTopSegment()->isSimpleInversion()) &&
     topSegment->getLength() <= _gapThreshold;*/
}


bool DefaultRearrangement::isGap(BottomSegmentIteratorConstPtr bottomSegment)
{
  return bottomSegment->hasChild(_childIndex) == false &&
     bottomSegment->getLength() <= _gapThreshold;
/*
  return 
     (bottomSegment->getBottomSegment()->isGapDeletion(_childIndex) ||
      bottomSegment->getBottomSegment()->isSimpleInversion(_childIndex)) &&
     bottomSegment->getLength() <= _gapThreshold;
*/
}

bool DefaultRearrangement::
isForwardAligned(TopSegmentIteratorConstPtr inLeft, 
                 TopSegmentIteratorConstPtr inRight,
                 BottomSegmentIteratorConstPtr inParentLeft,
                 BottomSegmentIteratorConstPtr inParentRight)
{
  TopSegmentIteratorConstPtr child = inLeft->copy();
  TopSegmentIteratorConstPtr temp = child->copy();
  BottomSegmentIteratorConstPtr parent = inParentLeft->copy();
  assert(child->getReversed() == false && parent->getReversed() == false);
  assert(inLeft->getTopSegment()->getArrayIndex() <=
         inRight->getTopSegment()->getArrayIndex());
  assert(inParentLeft->getBottomSegment()->getArrayIndex() <=
         inParentRight->getBottomSegment()->getArrayIndex());
  assert(isGap(inLeft) == false && isGap(inRight) == false);
  assert(isGap(inParentLeft) == false && isGap(inParentRight) == false);
  assert(inLeft->getTopSegment()->getSequence() == 
         inRight->getTopSegment()->getSequence());
  assert(inParentLeft->getBottomSegment()->getSequence() == 
         inParentRight->getBottomSegment()->getSequence());
  
  while (child->equals(inRight) == false)
  {
    temp->toChild(parent, _childIndex);
    if (child->equals(temp) == false || temp->getReversed() == true)
    {
      return false;
    }
    child->toRight();
    while (isGap(child))
    {
      child->toRight();
    }
    parent->toRight();
    while (isGap(parent))
    {
      parent->toRight();
    }    
    if (parent->equals(inParentRight) == true &&
        child->equals(inRight) == false)
    {
      return false;
    }
  }
  return parent->equals(inParentRight);
}

bool DefaultRearrangement::
isReverseAligned(TopSegmentIteratorConstPtr inLeft, 
                TopSegmentIteratorConstPtr inRight,
                BottomSegmentIteratorConstPtr inParentLeft,
                BottomSegmentIteratorConstPtr inParentRight)
{
  TopSegmentIteratorConstPtr child = inRight->copy();
  TopSegmentIteratorConstPtr temp = child->copy();
  BottomSegmentIteratorConstPtr parent = inParentLeft->copy();
  assert(child->getReversed() == false && parent->getReversed() == false);
  assert(inLeft->getTopSegment()->getArrayIndex() <=
         inRight->getTopSegment()->getArrayIndex());
  assert(inParentLeft->getBottomSegment()->getArrayIndex() <=
         inParentRight->getBottomSegment()->getArrayIndex());
  assert(isGap(inLeft) == false && isGap(inRight) == false);
  assert(isGap(inParentLeft) == false && isGap(inParentRight) == false);
  assert(inLeft->getTopSegment()->getSequence() == 
         inRight->getTopSegment()->getSequence());
  assert(inParentLeft->getBottomSegment()->getSequence() == 
         inParentRight->getBottomSegment()->getSequence());
  
  while (child->equals(inLeft) == false)
  {
    temp->toChild(parent, _childIndex);
    if (child->equals(temp) == false || temp->getReversed() == false)
    {
      return false;
    }
    child->toLeft();
    while (isGap(child))
    {
      child->toLeft();
    }
    parent->toRight();
    while (isGap(parent))
    {
      parent->toRight();
    }    
    if (parent->equals(inParentRight) == true && 
        child->equals(inLeft) == false)
    {
      return false;
    }
  }
  return parent->equals(inParentRight);
}

void DefaultRearrangement::findRight(TopSegmentIteratorConstPtr inLeft,
                                     TopSegmentIteratorConstPtr& outRight)
{
  outRight->copy(inLeft);
  if (outRight->getReversed() == true)
  {
    outRight->toReverse();
  }

  while (outRight->getTopSegment()->isLast() == false)
  {
    outRight->toRight();
    if (isGap(outRight) == false)
    {
      outRight->toLeft();
      break;
    }
    else
    {
      outRight->toRight();
    }
  }
  assert(outRight->getTopSegment()->getArrayIndex() >= 
         inLeft->getTopSegment()->getArrayIndex());
}

void DefaultRearrangement::findLeft(TopSegmentIteratorConstPtr inRight,
                                    TopSegmentIteratorConstPtr& outLeft)
{
  outLeft->copy(inRight);
  if (outLeft->getReversed() == true)
  {
    outLeft->toReverse();
  }

  while (outLeft->getTopSegment()->isFirst() == false)
  {
    outLeft->toLeft();
    if (isGap(outLeft) == false)
    {
      outLeft->toRight();
      break;
    }
    else
    {
      outLeft->toLeft();
    }
  }
  assert(inRight->getTopSegment()->getArrayIndex() >= 
         outLeft->getTopSegment()->getArrayIndex());
}

void DefaultRearrangement::findRight(BottomSegmentIteratorConstPtr inLeft,
                                     BottomSegmentIteratorConstPtr& outRight)
{
  outRight->copy(inLeft);
  if (outRight->getReversed() == true)
  {
    outRight->toReverse();
  }

  while (outRight->getBottomSegment()->isLast() == false)
  {
    outRight->toRight();
    if (isGap(outRight) == false)
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
                                    BottomSegmentIteratorConstPtr& outLeft)
{
  outLeft->copy(inRight);
  if (outLeft->getReversed() == true)
  {
    outLeft->toReverse();
  }

  while (outLeft->getBottomSegment()->isFirst() == false)
  {
    outLeft->toLeft();
    if (isGap(outLeft) == false)
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
                                  BottomSegmentIteratorConstPtr& outParentLeft,
                                  BottomSegmentIteratorConstPtr& outParentRight)
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
  assert(outParentLeft->getBottomSegment()->getArrayIndex() <=
         outParentRight->getBottomSegment()->getArrayIndex());


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
                                 TopSegmentIteratorConstPtr& outLeft,
                                 TopSegmentIteratorConstPtr& outRight)
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

// Segment is an inverted descendant of another Segment
// (Note that it can still be part of a more complex operation 
// such as a transposition -- which would have to be tested for
// in the insertion cycle code.  so this method is insufficient to
// return an ID=Inversion)
bool DefaultRearrangement::scanInversionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  resetStatus(topSegment);
  
  if (_startLeft->hasParent() == true && 
      _startLeft->getTopSegment()->getParentReversed() == true)
  {
    assert(_startRight->hasParent() == true);
    assert(_startRight->getTopSegment()->getParentReversed() == true);

    cout << "sl " << _startLeft << endl
         << "sr " << _startRight << endl
         << "pl " << _parentStartLeft << endl
         << "pr " << _parentStartRight << endl;

    findParents(_startLeft, _startRight, _parentStartLeft, _parentStartRight);

    cout << endl << endl;
    cout << "sl " << _startLeft << endl
         << "sr " << _startRight << endl
         << "pl " << _parentStartLeft << endl
         << "pr " << _parentStartRight << endl;

    _parentEndLeft->copy(_parentStartRight);
    _parentEndRight->copy(_parentStartRight);
    _startLeft2->toChild(_parentStartRight, _childIndex);
    _startRight2->toChild(_parentStartLeft, _childIndex);
    
    if (isReverseAligned(_startLeft, _startRight, 
                         _parentStartLeft, _parentStartRight) == true)
    {
      return true;
    }
  }

  return false;
}
