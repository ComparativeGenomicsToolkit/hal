/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hal.h"
#include "halTopSegmentIterator.h"
#include "defaultSegmentIterator.h"

using namespace std;
using namespace hal;

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE OVERRIDE
//////////////////////////////////////////////////////////////////////////////
void TopSegmentIterator::print(ostream& os) const
{
  os << "TopSegIt: ";
  DefaultSegmentIterator::print(os);

  hal_index_t ai = getArrayIndex();
  bool offRight = 
     isTop() ? ai >= (hal_index_t)getGenome()->getNumTopSegments() :
     ai >= (hal_index_t)getGenome()->getNumBottomSegments();

  if (ai != NULL_INDEX && !offRight)
  {
    os << " pIdx=" << getParentIndex() << " npIdx=" << getNextParalogyIndex();
  }
}

//////////////////////////////////////////////////////////////////////////////
// TOP SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
hal_index_t TopSegmentIterator::getParentIndex() const
{
  return _topSegment->getParentIndex();
}

bool TopSegmentIterator::hasParent() const
{
  assert(inRange() == true);
  return _topSegment->getParentIndex() != NULL_INDEX;
}

void TopSegmentIterator::setParentIndex(hal_index_t parIdx)
{
  _topSegment->setParentIndex(parIdx);
}

bool TopSegmentIterator::getParentReversed() const
{
  return _topSegment->getParentReversed();
}

void TopSegmentIterator::setParentReversed(bool isReversed)
{
  _topSegment->setParentReversed(isReversed);
}

hal_index_t TopSegmentIterator::getBottomParseIndex() const
{
  return _topSegment->getBottomParseIndex();
}

void TopSegmentIterator::setBottomParseIndex(hal_index_t botParseIdx)
{
  _topSegment->setBottomParseIndex(botParseIdx);
}

hal_offset_t TopSegmentIterator::getBottomParseOffset() const
{
  return _topSegment->getBottomParseOffset();
}

bool TopSegmentIterator::hasParseDown() const
{
  assert (inRange() == true);
  assert (_topSegment->getBottomParseIndex() == NULL_INDEX ||
          _topSegment->getGenome()->getNumChildren() > 0);
  return _topSegment->getBottomParseIndex() != NULL_INDEX;
}

hal_index_t TopSegmentIterator::getNextParalogyIndex() const
{
  return _topSegment->getNextParalogyIndex();
}

bool TopSegmentIterator::hasNextParalogy() const
{
  return _topSegment->hasNextParalogy();
}

void TopSegmentIterator::setNextParalogyIndex(hal_index_t parIdx)
{
  _topSegment->setNextParalogyIndex(parIdx);
}

hal_index_t TopSegmentIterator::getLeftParentIndex() const
{
  return _topSegment->getLeftParentIndex();
}

hal_index_t TopSegmentIterator::getRightParentIndex() const
{
  return _topSegment->getRightParentIndex();
}

bool TopSegmentIterator::isCanonicalParalog() const
{
  return _topSegment->isCanonicalParalog();
}

//////////////////////////////////////////////////////////////////////////////
// TOP SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void TopSegmentIterator::toChild(BottomSegmentIteratorConstPtr bs, 
                                     hal_size_t child) const
{
  _topSegment->setArrayIndex(bs->getGenome()->getChild(child),
                             bs->getChildIndex(child));
  _startOffset = bs->getStartOffset();
  _endOffset = bs->getEndOffset();
  _reversed = bs->getReversed();
  if (bs->getChildReversed(child) == true)
  {
    toReverse();
  }
  assert (inRange() == true);
}

void TopSegmentIterator::toChildG(BottomSegmentIteratorConstPtr bs, 
                                      const Genome* childGenome) const
{
  hal_index_t childIndex = bs->getGenome()->getChildIndex(childGenome);
  assert(childIndex != NULL_INDEX);
  toChild(bs, childIndex);
  assert (inRange() == true);
  assert(getGenome() == childGenome);
}


void TopSegmentIterator::toParseUp(BottomSegmentIteratorConstPtr bs) const
{ 
  const Genome* genome = bs->getGenome();
  hal_index_t index = bs->getTopParseIndex();

  _topSegment->setArrayIndex(genome, index);
  _reversed = bs->getReversed();

  hal_index_t startPos = bs->getStartPosition();

  while (startPos >= _topSegment->getStartPosition() + 
         (hal_index_t)_topSegment->getLength())
  {
    _topSegment->setArrayIndex(genome, ++index);
  }

  if (_reversed == false)
  {
    _startOffset = startPos - _topSegment->getStartPosition();    
    hal_index_t topEnd = _topSegment->getStartPosition() + 
       (hal_index_t)_topSegment->getLength();
    hal_index_t botEnd = bs->getStartPosition() + 
       (hal_index_t)bs->getLength();
    _endOffset = max((hal_index_t)0, topEnd - botEnd);
  }
  else
  {
    _startOffset = _topSegment->getStartPosition() + _topSegment->getLength() -
       1 - startPos;
    hal_index_t topEnd = _topSegment->getStartPosition();
    hal_index_t botEnd = bs->getStartPosition() - 
       (hal_index_t)bs->getLength() + 1;
    _endOffset = max((hal_index_t)0, botEnd - topEnd);  
  }
  assert (_startOffset + _endOffset <= _topSegment->getLength());
  assert (inRange() == true);
}

TopSegmentIteratorPtr TopSegmentIterator::copy()
{
  TopSegmentIteratorPtr newIt = 
     getGenome()->getTopSegmentIterator(getArrayIndex());
  if (_reversed)
  {
    newIt->toReverse();
  }
  newIt->slice(_startOffset, _endOffset);
  return newIt;
}

TopSegmentIteratorConstPtr TopSegmentIterator::copy() const
{
  TopSegmentIteratorConstPtr newIt = 
     getGenome()->getTopSegmentIterator(getArrayIndex());
  if (_reversed)
  {
    newIt->toReverse();
  }
  newIt->slice(_startOffset, _endOffset);
  return newIt;
}

void TopSegmentIterator::copy(TopSegmentIteratorConstPtr ts) const
{
  _topSegment->setArrayIndex(ts->getGenome(), ts->getArrayIndex());
  _startOffset = ts->getStartOffset();
  _endOffset = ts->getEndOffset();
  _reversed = ts->getReversed();
}

bool TopSegmentIterator::equals(TopSegmentIteratorConstPtr other) const
{
  assert(getGenome() == other->getGenome());
  return getArrayIndex() == other->getArrayIndex();
}

void TopSegmentIterator::toNextParalogy() const
{
  assert(_topSegment->getNextParalogyIndex() != NULL_INDEX);
  assert(_topSegment->getNextParalogyIndex() != _topSegment->getArrayIndex());
  bool rev = _topSegment->getParentReversed();
  _topSegment->setArrayIndex(getGenome(), _topSegment->getNextParalogyIndex());
  if(_topSegment->getParentReversed() != rev)
  {
    toReverse();
  }
}

