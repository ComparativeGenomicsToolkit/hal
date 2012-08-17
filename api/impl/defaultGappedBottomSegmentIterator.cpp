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
#include "hal.h"
#include "halGappedBottomSegmentIterator.h"
#include "defaultGappedBottomSegmentIterator.h"

using namespace std;
using namespace hal;

DefaultGappedBottomSegmentIterator::DefaultGappedBottomSegmentIterator(
  BottomSegmentIteratorConstPtr left,
  hal_size_t childIndex,
  hal_size_t gapThreshold) :
  _childIndex(childIndex),
  _gapThreshold(gapThreshold)
{
  if (left->getStartOffset() != 0 || left->getEndOffset() != 0)
  {
    throw hal_exception("offset not currently supported in gapped iterators");
  }
  const Genome* child = 
     left->getBottomSegment()->getGenome()->getChild(_childIndex);
  if (child == NULL)
  {
    throw hal_exception("can't init GappedBottomIterator with no child genome");
  }
  _left = left->copy();
  _right = left->copy();
  _temp = left->copy();
  _temp2 = left->copy();
  _leftChild = child->getTopSegmentIterator();
  _rightChild = _leftChild->copy();
  extendRight();
}

DefaultGappedBottomSegmentIterator::~DefaultGappedBottomSegmentIterator()
{

}

 // Gpped Segment Iterator methods
hal_size_t DefaultGappedBottomSegmentIterator::getGapThreshold() const
{
  return _gapThreshold;
}

hal_size_t DefaultGappedBottomSegmentIterator::getChildIndex() const
{
  return _childIndex;
}

hal_size_t DefaultGappedBottomSegmentIterator::getNumSegments() const
{
  return (hal_size_t)abs(_right->getBottomSegment()->getArrayIndex() - 
                         _left->getBottomSegment()->getArrayIndex() + 1);
}

hal_size_t DefaultGappedBottomSegmentIterator::getNumGapInsertions() const
{
  return 0;
}

hal_size_t DefaultGappedBottomSegmentIterator::getNumGapDeletions() const
{
  return 0;
}

hal_size_t DefaultGappedBottomSegmentIterator::getNumGapInsertedBases() const
{
  return 0;
}

hal_size_t DefaultGappedBottomSegmentIterator::getNumGapDeletedBases() const
{
  return 0;
}

bool DefaultGappedBottomSegmentIterator::isLast() const
{
  throw hal_exception("not imp");
  return false;
}

bool DefaultGappedBottomSegmentIterator::isFirst() const
{
  throw hal_exception("not imp");
  return false;
}

hal_index_t DefaultGappedBottomSegmentIterator::getLeftArrayIndex() const
{
  throw hal_exception("not imp");
}

hal_index_t DefaultGappedBottomSegmentIterator::getRightArrayIndex() const
{
  throw hal_exception("not imp");
}

// Segment Iterator methods
void DefaultGappedBottomSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));

  _right->copy(_left);
  _right->toLeft();
  _left->copy(_right);

  hal_index_t ns =
     _right->getBottomSegment()->getGenome()->getNumBottomSegments();
  hal_index_t i = _right->getBottomSegment()->getArrayIndex();
  if ((_right->getReversed() == true && i < ns) ||
      (_right->getReversed() == false && i > 0))
  {
    extendLeft();
  }
}

void DefaultGappedBottomSegmentIterator::toRight(hal_index_t rightCutoff) const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));

  _left->copy(_right);
  _left->toRight();
  _right->copy(_left);
  
  hal_index_t ns = 
     _left->getBottomSegment()->getGenome()->getNumBottomSegments();
  hal_index_t i = _left->getBottomSegment()->getArrayIndex();
  if ((_left->getReversed() == true && i > 0) ||
      (_left->getReversed() == false && i < ns))
  {
    extendRight();
  }
}

void DefaultGappedBottomSegmentIterator::toReverse() const
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->leftOf(_right->getStartPosition()));
  _left->toReverse();
  _right->toReverse();
}

void DefaultGappedBottomSegmentIterator::toSite(hal_index_t position, 
                                             bool slice) const
{
  throw hal_exception("tosite not currently supported in gapped iterators");
}

bool DefaultGappedBottomSegmentIterator::hasNextParalogy() const
{
  return false;
}

void DefaultGappedBottomSegmentIterator::toNextParalogy() const
{
  return;
}

hal_offset_t DefaultGappedBottomSegmentIterator::getStartOffset() const
{
  return 0;
}

hal_offset_t DefaultGappedBottomSegmentIterator::getEndOffset() const
{
  return 0;
}

void DefaultGappedBottomSegmentIterator::slice(hal_offset_t startOffset,
                                            hal_offset_t endOffset) const
{
  throw hal_exception("slice not currently supported in gapped iterators");
}

hal_index_t DefaultGappedBottomSegmentIterator::getStartPosition() const
{
  return 
     getReversed() ? _right->getStartPosition() : _left->getStartPosition();
}

hal_size_t DefaultGappedBottomSegmentIterator::getLength() const
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

hal_bool_t DefaultGappedBottomSegmentIterator::getReversed() const
{
  assert(_left->getReversed() == _right->getReversed());
  return _left->getReversed();
}

void DefaultGappedBottomSegmentIterator::getString(std::string& outString) const
{
  throw hal_exception("getString not currently supported in gapped iterators");
}

bool DefaultGappedBottomSegmentIterator::leftOf(hal_index_t genomePos) const
{
  return _right->leftOf(genomePos);
}

bool DefaultGappedBottomSegmentIterator::rightOf(hal_index_t genomePos) const
{
  return _left->rightOf(genomePos);
}

bool DefaultGappedBottomSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

// GappedBottomSegmentIterator Interface Methods
GappedBottomSegmentIteratorPtr DefaultGappedBottomSegmentIterator::copy()
{
  DefaultGappedBottomSegmentIterator* newIt =
     new DefaultGappedBottomSegmentIterator(_left, _childIndex, _gapThreshold);
  return GappedBottomSegmentIteratorPtr(newIt);
}

GappedBottomSegmentIteratorConstPtr DefaultGappedBottomSegmentIterator::copy() const
{
  const DefaultGappedBottomSegmentIterator* newIt =
     new DefaultGappedBottomSegmentIterator(_left, _childIndex, _gapThreshold);
  return GappedBottomSegmentIteratorConstPtr(newIt);
}

void DefaultGappedBottomSegmentIterator::copy(
  GappedBottomSegmentIteratorConstPtr ts) const
{
  _left->copy(ts->getLeft());
  _right->copy(ts->getRight());
  _childIndex = ts->getChildIndex();
  _gapThreshold = ts->getGapThreshold();
}

void DefaultGappedBottomSegmentIterator::toParent(
  GappedTopSegmentIteratorConstPtr ts) const
{
  TopSegmentIteratorConstPtr leftChild = ts->getLeft()->copy();
  TopSegmentIteratorConstPtr rightChild = ts->getRight()->copy();

  toRightNextUngapped(leftChild);
  toLeftNextUngapped(rightChild);

  _left->toParent(leftChild);
  _right->toParent(rightChild);
}

bool DefaultGappedBottomSegmentIterator::equals(
  GappedBottomSegmentIteratorConstPtr other) const
{
  return _left->equals(other->getLeft()) && _right->equals(other->getRight());
}

bool DefaultGappedBottomSegmentIterator::hasChild() const
{
  _temp->copy(_left);
  while(_temp != _right)
  {
    if (_temp->hasChild(_childIndex) == true)
    {
      return true;
    }
  }
  return false;
}

BottomSegmentIteratorConstPtr DefaultGappedBottomSegmentIterator::getLeft() const
{
  return _left;
}

BottomSegmentIteratorConstPtr DefaultGappedBottomSegmentIterator::getRight() const
{
  return _right;
}

// Internal methods

bool DefaultGappedBottomSegmentIterator::compatible(
  BottomSegmentIteratorConstPtr left,
  BottomSegmentIteratorConstPtr right) const
{
  assert(left->hasChild(_childIndex) && right->hasChild(_childIndex));
  assert(left->equals(right) == false);
  
  _leftChild->toChild(left, _childIndex);
  _rightChild->toChild(right, _childIndex);

  if (_leftChild->getTopSegment()->getParentReversed() != 
      _rightChild->getTopSegment()->getParentReversed())
  {
    return false;
  }

  if (_leftChild->leftOf(_leftChild->getStartPosition()) == false)
  {
    return false;
  }
  
  while (true)
  {
    assert(_leftChild->getTopSegment()->isLast() == false);
    _leftChild->toRight();
    if (_leftChild->hasParent() == true || 
        _leftChild->getLength() >= _gapThreshold)
    {
      if (_leftChild->equals(_rightChild))
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

void DefaultGappedBottomSegmentIterator::extendRight() const
{
  _right->copy(_left);
  if ((!_right->getReversed() && _right->getBottomSegment()->isLast()) ||
      (_right->getReversed() && _right->getBottomSegment()->isFirst()))
  {
    return;
  }
  _temp->copy(_right);

  while ((!_right->getReversed() && !_right->getBottomSegment()->isLast()) ||
         (_right->getReversed() && !_right->getBottomSegment()->isFirst()))
  {
    _right->toRight();
    if ((_right->hasChild(_childIndex) == false && 
         _right->getLength() > _gapThreshold) ||
        (_right->hasChild(_childIndex) == true && 
         compatible(_temp, _right) == false))
    {
      _right->toLeft();
      break;
    }
  }
}

void DefaultGappedBottomSegmentIterator::extendLeft() const
{
  _left->copy(_right);
  if ((!_left->getReversed() && _left->getBottomSegment()->isFirst()) ||
      (_left->getReversed() && _left->getBottomSegment()->isLast()))
  {
    return;
  }
  _temp->copy(_left);

  while ((!_left->getReversed() && !_left->getBottomSegment()->isLast()) ||
         (_left->getReversed() && !_left->getBottomSegment()->isFirst()))
  {
    _left->toLeft();
    if ((_left->hasChild(_childIndex) == false &&
         _left->getLength() > _gapThreshold) ||
        (_left->hasChild(_childIndex) == true &&
         compatible(_left, _temp) == false))
    {
      _left->toRight();
      break;
    }
  }
}

void DefaultGappedBottomSegmentIterator::toLeftNextUngapped(
  TopSegmentIteratorConstPtr ts) const
{
  while (ts->hasParent() == false && 
         ts->getLength() <= _gapThreshold)
  {
    if ((!ts->getReversed() && ts->getTopSegment()->isFirst() ||
         (ts->getReversed() && ts->getTopSegment()->isLast())))
    {
      break;
    }
    ts->toLeft();
  }
}

void DefaultGappedBottomSegmentIterator::toRightNextUngapped(
  TopSegmentIteratorConstPtr ts) const
{
  while (ts->hasParent() == false &&
         ts->getLength() <= _gapThreshold)
  {
    if ((!ts->getReversed() && ts->getTopSegment()->isLast() ||
         (ts->getReversed() && ts->getTopSegment()->isFirst())))
    {
      break;
    }
    ts->toRight();
  }
}
   
