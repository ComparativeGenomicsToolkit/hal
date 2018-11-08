/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <deque>
#include "halGenome.h"
#include "halDNAIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"

using namespace std;
using namespace hal;

GappedBottomSegmentIterator::GappedBottomSegmentIterator(
  BottomSegmentIteratorPtr left,
  hal_size_t childIndex,
  hal_size_t gapThreshold,
  bool atomic) :
  _childIndex(childIndex),
  _gapThreshold(gapThreshold),
  _atomic(atomic)
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
  assert(_atomic == false || _gapThreshold == 0);
  _left = left->clone();
  _right = left->clone();
  _temp = left->clone();
  _temp2 = left->clone();
  _leftChild = child->getTopSegmentIterator();
  _rightChild = _leftChild->clone();
  _leftDup = _leftChild->clone();
  _rightDup = _leftChild->clone();
  extendRight();
}

Segment* GappedBottomSegmentIterator::getSegment()
{
    return const_cast<Segment*>(_left->getSegment());

}

const Segment* GappedBottomSegmentIterator::getSegment() const
{
    return _left->getSegment();
}

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void GappedBottomSegmentIterator::setArrayIndex(Genome* genome, 
                                                hal_index_t arrayIndex)
{
  setLeft(genome->getBottomSegmentIterator(arrayIndex));
}

const Genome* GappedBottomSegmentIterator::getGenome() const
{
  return _left->getGenome();
}

Genome* GappedBottomSegmentIterator::getGenome()
{
  throw hal_exception("getGenome not supported in gapped iterators");
}

const Sequence* GappedBottomSegmentIterator::getSequence() const
{
  assert(_left->getBottomSegment()->getSequence() ==
         _right->getBottomSegment()->getSequence());
  return _left->getBottomSegment()->getSequence();
}

hal_index_t GappedBottomSegmentIterator::getStartPosition() const
{
  return _left->getStartPosition();
}

hal_index_t GappedBottomSegmentIterator::getEndPosition() const
{
  return _right->getEndPosition();
}

hal_size_t GappedBottomSegmentIterator::getLength() const
{
  return abs(getEndPosition() - getStartPosition()) + 1;
}

void GappedBottomSegmentIterator::getString(std::string& outString) const
{
  throw hal_exception("getString not supported in gapped iterators");
}

void GappedBottomSegmentIterator::setCoordinates(hal_index_t startPos, 
                                                        hal_size_t length)
{
  throw hal_exception("setCoordinates not supported in gapped iterators");
}

hal_index_t GappedBottomSegmentIterator::getArrayIndex() const
{
  throw hal_exception("getArrayIndex not supported in gapped iterators");
  return NULL_INDEX;
}

bool GappedBottomSegmentIterator::leftOf(hal_index_t genomePos) const
{
  return _right->leftOf(genomePos);
}

bool GappedBottomSegmentIterator::rightOf(hal_index_t genomePos) const
{
  return _left->rightOf(genomePos);
}

bool GappedBottomSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

bool GappedBottomSegmentIterator::isLast() const
{
  return _right->isLast();
}

bool GappedBottomSegmentIterator::isFirst() const
{
  return _left->isFirst();
}

bool GappedBottomSegmentIterator::isMissingData(double nThreshold) const
{
  if (nThreshold >= 1.0)
  {
    return false;
  }
  hal_index_t start = min(_left->getStartPosition(), _right->getEndPosition());
  DNAIteratorPtr di = _left->getGenome()->getDNAIterator(start);
  hal_size_t length = getLength();
  size_t maxNs = nThreshold * (double)length;
  size_t Ns = 0;
  char c;
  for (size_t i = 0; i < length; ++i, di->toRight())
  {
    c = di->getChar();
    if (c == 'N' || c == 'n')
    {
      ++Ns;
    }
    if (Ns > maxNs)
    {
      return true;
    }
    if ((length - i) < (maxNs - Ns))
    {
      break;
    }
  }
  return false;
}

bool GappedBottomSegmentIterator::isTop() const
{
  return false;
}

hal_size_t GappedBottomSegmentIterator::getMappedSegments(
  MappedSegmentSet& outSegments,
  const Genome* tgtGenome,
  const set<const Genome*>* genomesOnPath,
  bool doDupes,
  hal_size_t minLength,
  const Genome *coalescenceLimit,
  const Genome *mrca) const
{
  throw hal_exception("getMappedSegments is not supported in GappedTopSegmentIterator");
}

void GappedBottomSegmentIterator::print(std::ostream& os) const
{
  os << "Gapped Bottom Segment: (thresh=" << getGapThreshold() 
     << " ci=" << getChildIndex() << ")\n";
  os << "Left: ";
  if (_left.get() == NULL)
  {
    os << "NULL";
  }
  else
  {
    os << *_left;
  }
  os << "\nRight: ";
  if (_right.get() == NULL)
  {
    os << "NULL";
  }
  else
  {
    os << *_right;
  }
}

//////////////////////////////////////////////////////////////////////////////
// SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void GappedBottomSegmentIterator::toLeft(hal_index_t leftCutoff)
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->getReversed() || 
         _left->leftOf(_right->getStartPosition()));
  assert(_left->equals(_right) || !_left->getReversed() || 
         _left->rightOf(_right->getStartPosition()));

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

void GappedBottomSegmentIterator::toRight(hal_index_t rightCutoff)
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->getReversed() || 
         _left->leftOf(_right->getStartPosition()));
  assert(_left->equals(_right) || !_left->getReversed() || 
         _left->rightOf(_right->getStartPosition()));

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

void GappedBottomSegmentIterator::toReverse()
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->getReversed() || 
         _left->leftOf(_right->getStartPosition()));
  assert(_left->equals(_right) || !_left->getReversed() || 
         _left->rightOf(_right->getStartPosition()));

  _left->toReverse();
  _right->toReverse();
  swap(_left, _right);
}

void GappedBottomSegmentIterator::toReverseInPlace()
{
  assert(_right->getReversed() == _left->getReversed());
  assert(_left->equals(_right) || _left->getReversed() || 
         _left->leftOf(_right->getStartPosition()));
  assert(_left->equals(_right) || !_left->getReversed() || 
         _left->rightOf(_right->getStartPosition()));

  _left->toReverseInPlace();
  _right->toReverseInPlace();
  swap(_left, _right);
}

void GappedBottomSegmentIterator::toSite(hal_index_t position, 
                                             bool slice)
{
  throw hal_exception("tosite not currently supported in gapped iterators");
}

hal_offset_t GappedBottomSegmentIterator::getStartOffset() const
{
  return 0;
}

hal_offset_t GappedBottomSegmentIterator::getEndOffset() const
{
  return 0;
}

void GappedBottomSegmentIterator::slice(hal_offset_t startOffset,
                                            hal_offset_t endOffset)
{
  throw hal_exception("slice not currently supported in gapped iterators");
}

bool GappedBottomSegmentIterator::getReversed() const
{
  assert(_left->getReversed() == _right->getReversed());
  return _left->getReversed();
}

//////////////////////////////////////////////////////////////////////////////
// GAPPED SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
hal_size_t GappedBottomSegmentIterator::getGapThreshold() const
{
  return _gapThreshold;
}

bool GappedBottomSegmentIterator::getAtomic() const
{
  return _atomic;
}

hal_size_t GappedBottomSegmentIterator::getChildIndex() const
{
  return _childIndex;
}

hal_size_t GappedBottomSegmentIterator::getNumSegments() const
{
  return (hal_size_t)abs(_right->getBottomSegment()->getArrayIndex() - 
                         _left->getBottomSegment()->getArrayIndex()) + 1;
}

hal_size_t GappedBottomSegmentIterator::getNumGaps() const
{
  hal_size_t count = 0;
  _temp->copy(_left);
  for (; _temp->equals(_right) == false; _temp->toRight())
  {
    if (_temp->hasChild(_childIndex) == false && 
        _temp->getLength() <= _gapThreshold)
    {
      ++count;
    }
  }

  return count;
}

hal_size_t GappedBottomSegmentIterator::getNumGapBases() const
{
    hal_size_t count = 0;
  _temp->copy(_left);
  for (; _temp->equals(_right) == false; _temp->toRight())
  {
    if (_temp->hasChild(_childIndex) == false)
    {
      count += _temp->getLength();
    }
  }
  return count;
}

hal_index_t GappedBottomSegmentIterator::getLeftArrayIndex() const
{
  return _left->getBottomSegment()->getArrayIndex();
}

hal_index_t GappedBottomSegmentIterator::getRightArrayIndex() const
{
  return _right->getBottomSegment()->getArrayIndex();
}

//////////////////////////////////////////////////////////////////////////////
// GAPPED BOTTOM SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
GappedBottomSegmentIteratorPtr GappedBottomSegmentIterator::clone() const
{
  GappedBottomSegmentIterator* newIt =
     new GappedBottomSegmentIterator(_left, _childIndex, _gapThreshold,
       _atomic);
  newIt->_left->copy(_left);
  newIt->_right->copy(_right);
  newIt->_childIndex = _childIndex;
  newIt->_gapThreshold = _gapThreshold;
  assert(hasChild() == newIt->hasChild());
  assert(getChildReversed() == newIt->getChildReversed());
  return GappedBottomSegmentIteratorPtr(newIt);
}


void GappedBottomSegmentIterator::copy(
  GappedBottomSegmentIteratorPtr ts)
{
  _left->copy(ts->getLeft());
  _right->copy(ts->getRight());
  _childIndex = ts->getChildIndex();
  _gapThreshold = ts->getGapThreshold();
  assert(hasChild() == ts->hasChild());
  assert(getChildReversed() == ts->getChildReversed());
}

void GappedBottomSegmentIterator::toParent(
  GappedTopSegmentIteratorPtr ts)
{
  TopSegmentIteratorPtr leftChild = ts->getLeft()->clone();
  TopSegmentIteratorPtr rightChild = ts->getRight()->clone();

  const Genome* child = ts->getLeft()->getTopSegment()->getGenome();
  const Genome* parent = child->getParent();
  _childIndex = parent->getChildIndex(child);

  toRightNextUngapped(leftChild);
  toLeftNextUngapped(rightChild);

  _left->toParent(leftChild);
  _right->toParent(rightChild);

  // should extend here?!
}

bool GappedBottomSegmentIterator::equals(
  GappedBottomSegmentIteratorPtr other) const
{
  _temp->copy(_left);
  toRightNextUngapped(_temp);
  _temp2->copy(other->getLeft());
  toRightNextUngapped(_temp2);
  if (_temp->equals(_temp2) == false)
  {
    return false;
  }
  _temp->copy(_right);
  toLeftNextUngapped(_temp);
  _temp2->copy(other->getRight());
  toLeftNextUngapped(_temp2);
  return _temp->equals(_temp2);
}

bool GappedBottomSegmentIterator::adjacentTo(
  GappedBottomSegmentIteratorPtr other) const
{
  _temp->copy(_left);
  toLeftNextUngapped(_temp);
  if (_temp->isFirst() == false)
  {
    _temp2->copy(other->getLeft());
    if (_temp2->isFirst() == false)
    {
      _temp2->toLeft();
      toLeftNextUngapped(_temp2);
      if (_temp->equals(_temp2))
      {
        return true;
      }
    }
    _temp2->copy(other->getRight());
    if (_temp2->isLast() == false)
    {
      _temp2->toRight();
      toRightNextUngapped(_temp2);
      if (_temp->equals(_temp2))
      {
        return true;
      }
    }
  }

  _temp->copy(_right);
  toRightNextUngapped(_temp);
  if (_temp->isLast() == false)
  {
    _temp2->copy(other->getLeft());
    if (_temp2->isFirst() == false)
    {
      _temp2->toLeft();
      toLeftNextUngapped(_temp2);
      if (_temp->equals(_temp2))
      {
        return true;
      }
    }
    _temp2->copy(other->getRight());
    if (_temp2->isLast() == false)
    {
      _temp2->toRight();
      toRightNextUngapped(_temp2);
      if (_temp->equals(_temp2))
      {
        return true;
      }
    }
  }

  return false;
}

bool GappedBottomSegmentIterator::hasChild() const
{
  _temp->copy(_left);
  _temp2->copy(_right);
  
  toRightNextUngapped(_temp);
  toLeftNextUngapped(_temp2);

  if (_left->equals(_right))
  {
    //assert(_temp->equals(_left) && _temp2->equals(_right));
  }
  // to do: verify edge cases
  // assert(_temp->hasChild(_childIndex) == _temp2->hasChild(_childIndex));
  return _temp->hasChild(_childIndex) && _temp2->hasChild(_childIndex);
}

bool GappedBottomSegmentIterator::getChildReversed() const
{
  if (hasChild())
  {
    _temp->copy(_left);
    _temp2->copy(_right);

    toRightNextUngapped(_temp);
    toLeftNextUngapped(_temp2);
    
    if (_left->equals(_right))
    {
      assert(_temp->equals(_left) && _temp2->equals(_right));
    }

    assert(_temp->hasChild(_childIndex) && _temp2->hasChild(_childIndex));
    assert(_temp->getBottomSegment()->getChildReversed(_childIndex) == 
           _temp2->getBottomSegment()->getChildReversed(_childIndex));
    return _temp->getBottomSegment()->getChildReversed(_childIndex);
  }
  return false;
}

BottomSegmentIteratorPtr GappedBottomSegmentIterator::getLeft() const
{
  return _left;
}

BottomSegmentIteratorPtr GappedBottomSegmentIterator::getRight() const
{
  return _right;
}

void 
GappedBottomSegmentIterator::setLeft(BottomSegmentIteratorPtr bi)
{
  throw hal_exception("setLeft not currently supported in bottom iterator");
}

//////////////////////////////////////////////////////////////////////////////
// INTERNAL METHODS
//////////////////////////////////////////////////////////////////////////////
bool GappedBottomSegmentIterator::compatible(
  BottomSegmentIteratorPtr left,
  BottomSegmentIteratorPtr right) const
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

  if (_leftChild->hasNextParalogy() != _rightChild->hasNextParalogy())
  {
    return false;
  }
  if ((!_leftChild->getReversed() && 
       _leftChild->leftOf(_rightChild->getStartPosition()) == false) ||
      (_leftChild->getReversed() && 
       _leftChild->rightOf(_rightChild->getStartPosition()) == false))
  {    
    return false;
  }
  
  if (left->getBottomSegment()->getSequence() != 
      right->getBottomSegment()->getSequence() ||
      _leftChild->getTopSegment()->getSequence() != 
      _rightChild->getTopSegment()->getSequence())
  {
    return false;
  }
  
  while (true)
  {
    assert(_leftChild->isLast() == false);
    _leftChild->toRight();
    if (_leftChild->hasParent() == true || 
        _leftChild->getLength() > _gapThreshold)
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

  _leftChild->toChild(left, _childIndex);
  _rightChild->toChild(right, _childIndex);
  if (_leftChild->hasNextParalogy() == true)
  {
    _leftDup->copy(_leftChild);
    _leftDup->toNextParalogy();
    _rightDup->copy(_rightChild);
    _rightDup->toNextParalogy();
  
    if ((_leftDup->getReversed() == false && 
         _leftDup->leftOf(_rightDup->getStartPosition()) == false) ||
        (_leftDup->getReversed() == true && 
         _rightDup->leftOf(_leftDup->getStartPosition()) == false))
    {
      return false;
    }
    if (_leftDup->getTopSegment()->getSequence() != 
        _rightDup->getTopSegment()->getSequence())
    {
      return false;
    }

    while (true)
    {
      assert(_leftDup->isLast() == false);
      _leftDup->toRight();
      if (_leftDup->hasParent() == true || 
          _leftDup->getLength() > _gapThreshold)
      {
        if (_leftDup->equals(_rightDup))
        {
          break;
        }
        else
        {
          return false;
        }
      }
    }
  }

  return true;
}

void GappedBottomSegmentIterator::extendRight()
{
  _right->copy(_left);
  if (_atomic == true ||
      (!_right->getReversed() && _right->getBottomSegment()->isLast()) ||
      (_right->getReversed() && _right->getBottomSegment()->isFirst()))
  {
    return;
  }
  toRightNextUngapped(_right);
  _temp->copy(_right);

  while (_right->isLast() == false)
  {
    _right->toRight();
    toRightNextUngapped(_right);    
  
    if ((_right->hasChild(_childIndex) == false && 
         _right->getLength() > _gapThreshold) ||
        (_temp->hasChild(_childIndex) == false && 
         _temp->getLength() > _gapThreshold) ||
        (_right->hasChild(_childIndex) == true &&
         _temp->hasChild(_childIndex) == true &&
         compatible(_temp, _right) == false))
    {
      _right->toLeft();
      break;
    }
    _temp->toRight();
    toRightNextUngapped(_temp);
  }
}

void GappedBottomSegmentIterator::extendLeft()
{
  _left->copy(_right);
  if (_atomic == true ||
      (!_left->getReversed() && _left->getBottomSegment()->isFirst()) ||
      (_left->getReversed() && _left->getBottomSegment()->isLast()))
  {
    return;
  }
  toLeftNextUngapped(_left);
  _temp->copy(_left);
  
  while (_left->isFirst() == false)
  {
    _left->toLeft();
    toLeftNextUngapped(_left);

    if ((_left->hasChild(_childIndex) == false &&
         _left->getLength() > _gapThreshold) ||
        (_temp->hasChild(_childIndex) == false &&
         _temp->getLength() > _gapThreshold) ||
        (_left->hasChild(_childIndex) == true &&
         _temp->hasChild(_childIndex) == true && 
         compatible(_left, _temp) == false))
    {
      _left->toRight();
      break;
    }
    _temp->toLeft();
    toLeftNextUngapped(_temp);
  }
}

void GappedBottomSegmentIterator::toLeftNextUngapped(
  BottomSegmentIteratorPtr bs) const
{
    // FIXME: these should operate on self and copies be made
  while (bs->hasChild(_childIndex) == false && 
         bs->getLength() <= _gapThreshold)
  {
    if ((!bs->getReversed() && bs->getBottomSegment()->isFirst()) ||
         (bs->getReversed() && bs->getBottomSegment()->isLast()))
    {
      break;
    }
    bs->toLeft();
  }
}

void GappedBottomSegmentIterator::toRightNextUngapped(
  BottomSegmentIteratorPtr bs) const
{
  while (bs->hasChild(_childIndex) == false &&
         bs->getLength() <= _gapThreshold)
  {
    if ((!bs->getReversed() && bs->getBottomSegment()->isLast()) ||
         (bs->getReversed() && bs->getBottomSegment()->isFirst()))
    {
      break;
    }
    bs->toRight();
  }
}

void GappedBottomSegmentIterator::toLeftNextUngapped(
  TopSegmentIteratorPtr ts) const
{
  while (ts->hasParent() == false && 
         ts->getLength() <= _gapThreshold)
  {
    if ((!ts->getReversed() && ts->getTopSegment()->isFirst()) ||
         (ts->getReversed() && ts->getTopSegment()->isLast()))
    {
      break;
    }
    ts->toLeft();
  }
}

void GappedBottomSegmentIterator::toRightNextUngapped(
  TopSegmentIteratorPtr ts) const
{
  while (ts->hasParent() == false &&
         ts->getLength() <= _gapThreshold)
  {
    if ((!ts->getReversed() && ts->getTopSegment()->isLast()) ||
         (ts->getReversed() && ts->getTopSegment()->isFirst()))
    {
      break;
    }
    ts->toRight();
  }
}
   
