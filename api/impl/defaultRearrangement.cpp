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
  _id(Invalid)
{
  _parent = childGenome->getParent();
  assert(_parent != NULL);
  // just allocating here -- need to be properly init
  _cur = _genome->getGappedTopSegmentIterator(0, _gapThreshold);
  _next = _cur->copy();
  _left = _cur->copy();
  _right = _left->copy();
  _leftParent = _parent->getGappedBottomSegmentIterator(0, 0, _gapThreshold);
  _rightParent = _leftParent->copy();
  _curParent = _leftParent->copy();
  _nextParent = _leftParent->copy();
  _top = _cur->getLeft()->copy();
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
  return _id == Deletion ? _leftParent->getLength() : _cur->getLength();
}

hal_size_t DefaultRearrangement::getNumContainedGaps() const
{
  return _id == Deletion ? _leftParent->getNumGaps() : _cur->getNumGaps();
}

hal_size_t DefaultRearrangement::getNumContainedGapBases() const
{
  return
     _id == Deletion ? _leftParent->getNumGapBases() : _cur->getNumGapBases();
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
  if (scanNothingCycle(topSegment) == true)
  {
    _id = Nothing;
  }
  else if (scanInversionCycle(topSegment) == true)
  {
    _id = Inversion;
  }
  else if (scanInsertionCycle(topSegment) == true)
  {
    if (_cur->getLength() > _gapThreshold)
    {
      _id = _cur->hasParent() ? Transposition : Insertion;
    }
    else
    {
      _id = Gap;
    }
  }
  else if (scanDeletionCycle(topSegment) == true) 
  {
    if (_leftParent->getLength() > _gapThreshold &&
        _leftParent->hasChild() == false)
    {
      _id = Deletion;
    }           
    else if (_leftParent->getLength() <= _gapThreshold)
    {
      _id = Gap;
    }
    else
    {
      _id = Complex;
    }
  }
  else if (scanDuplicationCycle(topSegment) == true)
  {
    _id = Duplication;
  }
  else
  {
    resetStatus(topSegment);
    if (_cur->isFirst() == false && _cur->isLast() == false)
    {
      _left->toLeft();
      _right->toRight();
      if (_left->hasParent() == true && _right->hasParent() == true && 
          _cur->hasParent() == true)
      {
        _curParent->toParent(_cur);
        _leftParent->toParent(_left);
        _rightParent->toParent(_right);
/*
        cout << "\n\ncomplex\n"
             << "cur " << _cur 
             << "\nleft ** " << _left
             << "\nright ** " << _right
             << "\ntp ** " << _curParent
             << "\nlp ** " << _leftParent 
             << "\nrp ** " << _rightParent << endl;
*/
/*
    if (_cur->getLeft()->getTopSegment()->getArrayIndex() == 1219)
    {
      GappedBottomSegmentIteratorConstPtr bt = _leftParent->copy();
      cout << "bt copy " << bt << endl;
      bt->toRight();
      cout << "bt right " << bt << endl;
      _leftParent->toRight();
      cout << "\n after lp to right ** " << _leftParent;
      _rightParent->toLeft();
      cout << "\n after rp to left ** " << _rightParent << endl;
      scanDeletionCycle(topSegment);
//      exit(10);
    }
*/
      }
    }

    _id = Complex;
  }
  
  return true;
}

bool DefaultRearrangement::identifyNext()
{
  assert(_cur->getReversed() == false);
  // don't like this.  need to refactor interface to make better
  // use of gapped iterators
  _top->copy(_cur->getRight());
  if (_top->getTopSegment()->getArrayIndex() < _genome->getNumTopSegments() - 1)
  {
    _top->toRight();
    identifyFromLeftBreakpoint(_top);
    return true;
  }
  else
  {
    return false;
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
  assert(topSegment.get());
  _genome = topSegment->getTopSegment()->getGenome();
  _parent = _genome->getParent();
  assert(_parent != NULL);

  _cur->setLeft(topSegment);
  _next->copy(_cur);
  _left->copy(_cur);
  _right->copy(_left);
}

// Segment corresponds to no rearrangemnt.  This will happen when 
// there is a rearrangement in the homolgous segment in its sibling 
// genome.  In general, we can expect about half of segments to correspond
// to such cases.   
bool DefaultRearrangement::scanNothingCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  bool first = _cur->isFirst();
  bool last = _cur->isLast();

  if (_cur->hasParent() == false)
  {
    return false;
  }
  _curParent->toParent(_cur);
  if (first == false)
  {
    _left->toLeft();
    if (_left->hasParent() == false)
    {
      return false;
    }
    _leftParent->toParent(_left);
    if (_leftParent->adjacentTo(_curParent) == false)
    {
      return false;
    }
    if (_left->getParentReversed() == true)
    {
      if (_cur->getParentReversed() == false || 
          _leftParent->rightOf(_curParent->getStartPosition()) == false)
      {
         return false;
      }
    }
    else
    {
      if (_cur->getParentReversed() == true || 
          _leftParent->leftOf(_curParent->getStartPosition()) == false)
      {
         return false;
      }
    }
  }
  if (last == false)
  {
    _right->toRight();
    if (_right->hasParent() == false)
    {
      return false;
    }
    _rightParent->toParent(_right);
    if (_rightParent->adjacentTo(_curParent) == false)
    {
      return false;
    }
    if (_right->getParentReversed() == true)
    {
      if (_cur->getParentReversed() == false || 
          _rightParent->leftOf(_curParent->getStartPosition()) == false)
      {
         return false;
      }
    }
    else
    {
      if (_cur->getParentReversed() == true || 
          _rightParent->rightOf(_curParent->getStartPosition()) == false)
      {
         return false;
      }
    }
  }
  return last && first ? _cur->getParentReversed() : true;
}

// Segment is an inverted descendant of another Segment but 
// otherwise no rearrangement.  
bool DefaultRearrangement::scanInversionCycle(
  TopSegmentIteratorConstPtr topSegment)
{
  assert(topSegment.get());
  resetStatus(topSegment);
  bool first = _cur->isFirst();
  bool last = _cur->isLast();

  if (_cur->hasParent() == false)
  {
    return false;
  }
  _curParent->toParent(_cur);
  if (first == false)
  {
    _left->toLeft();
    if (_left->hasParent() == false)
    {
      return false;
    }
    _leftParent->toParent(_left);
    if (_leftParent->adjacentTo(_curParent) == false)
    {
      return false;
    }
  }
  if (last == false)
  {
    _right->toRight();
    if (_right->hasParent() == false)
    {
      return false;
    }
    _rightParent->toParent(_right);
    if (_rightParent->adjacentTo(_curParent) == false)
    {
      return false;
    }
  }
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

  // eat up any adjacent insertions so they don't get double counted
  while (_next->hasParent() == false && _next->isLast() == false)
  {
    _right->copy(_next);
    _right->toRight();
    if (_right->hasParent() == false)
    {
      _next->copy(_right);
    }
    else
    {
      break;
    }
  }
  _right->copy(_next);
  assert(_next->equals(_cur) || _next->hasParent() == false);

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
    if (_left->hasParent() == true && _right->hasParent() == true)
    {
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
    _curParent->toParent(_right);
    return _curParent->equals(_rightParent);
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
