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

#include "defaultGappedTopSegmentIterator.h"
#include "hal.h"
using namespace std;
using namespace hal;

DefaultGappedTopSegmentIterator::DefaultGappedTopSegmentIterator(
  TopSegmentIteratorConstPtr left,
  hal_size_t childIndex,
  hal_size_t gapThreshold) :
  _childIndex(childIndex),
  _gapThreshold(gapThreshold)
{
  if (left->getStartOffset() != 0 || left->getEndOffset() != 0)
  {
    throw hal_exception("offset not currently supported in gapped iterators");
  }
  const Genome* parent = left->getTopSegment()->getGenome()->getParent();
  if (parent == NULL)
  {
    throw hal_exception("can't init GappedTopIterator with no parent genome");
  }
  _left = left->copy();
  _right = left->copy();
  _temp = left->copy();
  _temp2 = left->copy();
  _leftParent = parent->getBottomSegmentIterator();
  _rightParent = _leftParent->copy();
  extendRight();
}

DefaultGappedTopSegmentIterator::~DefaultGappedTopSegmentIterator()
{

}

 // Gpped Segment Iterator methods
hal_size_t DefaultGappedTopSegmentIterator::getGapThreshold() const
{
  return _gapThreshold;
}

hal_size_t DefaultGappedTopSegmentIterator::getNumSegments() const
{
  return (hal_size_t)abs(_right->getTopSegment()->getArrayIndex() - 
                         _left->getTopSegment()->getArrayIndex() + 1);
}

hal_size_t DefaultGappedTopSegmentIterator::getNumGapInsertions() const
{
  return 0;
}

hal_size_t DefaultGappedTopSegmentIterator::getNumGapDeletions() const
{
  return 0;
}

hal_size_t DefaultGappedTopSegmentIterator::getNumGapInsertedBases() const
{
  return 0;
}

hal_size_t DefaultGappedTopSegmentIterator::getNumGapDeletedBases() const
{
  return 0;
}

bool DefaultGappedTopSegmentIterator::isLast() const
{
  throw hal_exception("not imp");
  return false;
}

bool DefaultGappedTopSegmentIterator::isFirst() const
{
  throw hal_exception("not imp");
  return false;
}

hal_index_t DefaultGappedTopSegmentIterator::getLeftArrayIndex() const
{
  throw hal_exception("not imp");
}

hal_index_t DefaultGappedTopSegmentIterator::getRightArrayIndex() const
{
  throw hal_exception("not imp");
}

// Segment Iterator methods
void DefaultGappedTopSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));

  _right->copy(_left);
  _right->toLeft();
  extendLeft();
}

void DefaultGappedTopSegmentIterator::toRight(hal_index_t rightCutoff) const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));

  _left->copy(_right);
  _left->toRight();
  extendRight();
}

void DefaultGappedTopSegmentIterator::toReverse() const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));
  _left->toReverse();
  _right->toReverse();
}

void DefaultGappedTopSegmentIterator::toSite(hal_index_t position, 
                                             bool slice) const
{
  throw hal_exception("tosite not currently supported in gapped iterators");
}

bool DefaultGappedTopSegmentIterator::hasNextParalogy() const
{
  _temp->copy(_left);
  while(_temp != _right)
  {
    if (_temp->hasNextParalogy() == true || _temp->hasParent() == true)
    {
      return true;
    }
  }
  return false;
}

void DefaultGappedTopSegmentIterator::toNextParalogy() const
{
  _left->toNextParalogy();
  _right->toNextParalogy();
}

hal_offset_t DefaultGappedTopSegmentIterator::getStartOffset() const
{
  return 0;
}

hal_offset_t DefaultGappedTopSegmentIterator::getEndOffset() const
{
  return 0;
}

void DefaultGappedTopSegmentIterator::slice(hal_offset_t startOffset,
                                            hal_offset_t endOffset) const
{
  throw hal_exception("slice not currently supported in gapped iterators");
}

hal_index_t DefaultGappedTopSegmentIterator::getStartPosition() const
{
  return 
     getReversed() ? _right->getStartPosition() : _left->getStartPosition();
}

hal_size_t DefaultGappedTopSegmentIterator::getLength() const
{
  if (getReversed() == false)
  {
    return _right->getStartPosition() + _right->getLength() - 
       _left->getStartPosition();
  }
  else
  {
    return _left->getStartPosition() - _right->getStartPosition() + 
       _right->getLength();
  }
}

hal_bool_t DefaultGappedTopSegmentIterator::getReversed() const
{
  assert(_left->getReversed() == _right->getReversed());
  return _left->getReversed();
}

void DefaultGappedTopSegmentIterator::getString(std::string& outString) const
{
  throw hal_exception("getString not currently supported in gapped iterators");
}

bool DefaultGappedTopSegmentIterator::leftOf(hal_index_t genomePos) const
{
  return _right->leftOf(genomePos);
}

bool DefaultGappedTopSegmentIterator::rightOf(hal_index_t genomePos) const
{
  return _left->rightOf(genomePos);
}

bool DefaultGappedTopSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

// GappedTopSegmentIterator Interface Methods
GappedTopSegmentIteratorPtr DefaultGappedTopSegmentIterator::copy()
{
  DefaultGappedTopSegmentIterator* newIt =
     new DefaultGappedTopSegmentIterator(_left, _childIndex, _gapThreshold);
  return GappedTopSegmentIteratorPtr(newIt);
}

GappedTopSegmentIteratorConstPtr DefaultGappedTopSegmentIterator::copy() const
{
  const DefaultGappedTopSegmentIterator* newIt =
     new DefaultGappedTopSegmentIterator(_left, _childIndex, _gapThreshold);
  return GappedTopSegmentIteratorConstPtr(newIt);
}

void DefaultGappedTopSegmentIterator::copy(
  GappedTopSegmentIteratorConstPtr ts) const
{
  _left->copy(ts->getLeft());
  _right->copy(ts->getRight());
  _childIndex = ts->getChildIndex();
  _gapThreshold = ts->getGapThreshold();
}

void DefaultGappedTopSegmentIterator::toChild(
  GappedBottomSegmentIteratorConstPtr bs) const
{
  BottomSegmentIteratorConstPtr leftParent = bs->getLeft()->copy();
  BottomSegmentIteratorConstPtr rightParent = bs->getRight()->copy();

  toRightNextUngapped(leftParent);
  toLeftNextUngapped(rightParent);

  _left->toChild(leftParent, _childIndex);
  _right->toChild(rightParent, _childIndex);
}

bool DefaultGappedTopSegmentIterator::equals(
  GappedTopSegmentIteratorConstPtr other) const
{
  return _left->equals(other->getLeft()) && _right->equals(other->getRight());
}

bool DefaultGappedTopSegmentIterator::hasParent() const
{
  _temp->copy(_left);
  while(_temp != _right)
  {
    if (_temp->hasParent() == true)
    {
      return true;
    }
  }
  return false;
}

TopSegmentIteratorConstPtr DefaultGappedTopSegmentIterator::getLeft() const
{
  return _left;
}

TopSegmentIteratorConstPtr DefaultGappedTopSegmentIterator::getRight() const
{
  return _right;
}

// Internal methods

bool DefaultGappedTopSegmentIterator::compatible(
  TopSegmentIteratorConstPtr left,
  TopSegmentIteratorConstPtr right) const
{
  assert(left->hasParent() && right->hasParent());
  assert(left->equals(right) == false);
  
  _leftParent->toParent(left);
  _rightParent->toParent(right);

  if (_leftParent->getReversed() != _rightParent->getReversed())
  {
    return false;
  }

  if (_left->hasNextParalogy() != _right->hasNextParalogy())
  {
    return false;
  }
  if (_left->hasNextParalogy() == true && 
      _left->getTopSegment()->getNextParalogyReversed() !=
      _right->getTopSegment()->getNextParalogyReversed())
  {
    return false;
  }

  if (_leftParent->leftOf(_rightParent->getStartPosition()) == false)
  {
    return false;
  }
  
  while (true)
  {
    assert(_leftParent->getBottomSegment()->isLast() == false);
    _leftParent->toRight();
    if (_leftParent->hasChild(_childIndex) == true || 
        _leftParent->getLength() >= _gapThreshold)
    {
      if (_leftParent->equals(_rightParent))
      {
        break;
      }
      else
      {
        return false;
      }
    }
  }
  _temp->copy(_left);
  _temp->toNextParalogy();
  _temp2->copy(_right);
  _temp2->toNextParalogy();
  
  if (_temp->leftOf(_temp2->getStartPosition()) == false)
  {
    return false;
  }
  
  while (true)
  {
    assert(_temp->getTopSegment()->isLast() == false);
    _temp->toRight();
    if (_temp->hasParent() == true || 
        _temp->getLength() >= _gapThreshold)
    {
      if (_temp->equals(_temp2))
      {
        break;
      }
      else
      {
        return false;
      }
    }
  }

  return true;
}

void DefaultGappedTopSegmentIterator::extendRight() const
{
  _right->copy(_left);
  if ((!_left->getReversed() && _left->getTopSegment()->isLast()) ||
      (_left->getReversed() && _left->getTopSegment()->isFirst()))
  {
    return;
  }
  _right->toRight();
  _temp->copy(_right);

  while ((!_left->getReversed() && !_left->getTopSegment()->isLast()) ||
         (_left->getReversed() && !_left->getTopSegment()->isFirst()))
  {
    _right->toRight();
    if (_right->hasParent() == true || _right->getLength() > _gapThreshold)
    {
      if (compatible(_temp, _right) == false)
      {
        _right->toLeft();
        break;
      }
    }
  }
}

void DefaultGappedTopSegmentIterator::extendLeft() const
{
  _left->copy(_right);
  if ((!_right->getReversed() && _right->getTopSegment()->isFirst()) ||
      (_right->getReversed() && _right->getTopSegment()->isLast()))
  {
    return;
  }
  _left->toLeft();
  _temp->copy(_left);

  while ((!_left->getReversed() && !_left->getTopSegment()->isLast()) ||
         (_left->getReversed() && !_left->getTopSegment()->isFirst()))
  {
    _left->toLeft();
    if (_left->hasParent() == true || _left->getLength() > _gapThreshold)
    {
      if (compatible(_left, _temp) == false)
      {
        _left->toRight();
        break;
      }
    }
  }
}

void DefaultGappedTopSegmentIterator::toLeftNextUngapped(
  BottomSegmentIteratorConstPtr bs) const
{
  while (bs->hasChild(_childIndex) == false && 
         bs->getLength() <= _gapThreshold)
  {
    if ((!bs->getReversed() && bs->getBottomSegment()->isFirst() ||
         (bs->getReversed() && bs->getBottomSegment()->isLast())))
    {
      break;
    }
    bs->toLeft();
  }
}

void DefaultGappedTopSegmentIterator::toRightNextUngapped(
  BottomSegmentIteratorConstPtr bs) const
{
  while (bs->hasChild(_childIndex) == false &&
         bs->getLength() <= _gapThreshold)
  {
    if ((!bs->getReversed() && bs->getBottomSegment()->isLast() ||
         (bs->getReversed() && bs->getBottomSegment()->isFirst())))
    {
      break;
    }
    bs->toRight();
  }
}
   
