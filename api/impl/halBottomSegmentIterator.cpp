/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "halBottomSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include "halGenome.h"

using namespace std;
using namespace hal;

BottomSegmentIterator::
BottomSegmentIterator(BottomSegment* bottomSegment, 
                      hal_size_t startOffset, 
                      hal_size_t endOffset,
                      bool reversed) :
    SegmentIterator(startOffset, endOffset, reversed),
  _bottomSegment(bottomSegment)
{

}

hal_size_t BottomSegmentIterator::getNumSegmentsInGenome() const
{
  return getGenome()->getNumBottomSegments();
}
 
//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE OVERRIDE
//////////////////////////////////////////////////////////////////////////////
void BottomSegmentIterator::print(ostream& os) const
{
  os << "BotSegIt: ";
  SegmentIterator::print(os);

  hal_index_t ai = getArrayIndex();
  bool offRight = 
     isTop() ? ai >= (hal_index_t)getGenome()->getNumTopSegments() :
     ai >= (hal_index_t)getGenome()->getNumBottomSegments();
  
  if (ai != NULL_INDEX && !offRight)
  {
    os << " numChilds=" << getNumChildren();
    for (hal_size_t i = 0; i < getNumChildren(); ++i)
    {
      os << " cI[" << i << "]=" << getChildIndex(i);
      os << " cR[" << i << "]=" << getChildReversed(i);
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////
// BOTTOM SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
hal_size_t BottomSegmentIterator::getNumChildren() const
{
  return _bottomSegment->getNumChildren();
}

hal_index_t BottomSegmentIterator::getChildIndex(hal_size_t i) const
{
  return _bottomSegment->getChildIndex(i);
}

hal_index_t BottomSegmentIterator::getChildIndexG(
  const Genome* childGenome) const
{
  return _bottomSegment->getChildIndexG(childGenome);
}

bool BottomSegmentIterator::hasChild(hal_size_t child) const
{
  return _bottomSegment->hasChild(child);
}

bool BottomSegmentIterator::hasChildG(const Genome* childGenome) const
{
  return _bottomSegment->hasChildG(childGenome);
}

void BottomSegmentIterator::setChildIndex(hal_size_t i, 
                                              hal_index_t childIndex)
{
  _bottomSegment->setChildIndex(i, childIndex);
}

bool BottomSegmentIterator::getChildReversed(hal_size_t i) const
{
  return _bottomSegment->getChildReversed(i);
}

void BottomSegmentIterator::setChildReversed(hal_size_t child,
                                                 bool isReversed)
{
  _bottomSegment->setChildReversed(child, isReversed);
}

hal_index_t BottomSegmentIterator::getTopParseIndex() const
{
  return _bottomSegment->getTopParseIndex();
}

void BottomSegmentIterator::setTopParseIndex(hal_index_t parseIndex)
{
  _bottomSegment->setTopParseIndex(parseIndex);
}

hal_offset_t BottomSegmentIterator::getTopParseOffset() const
{
  return _bottomSegment->getTopParseOffset();
}

bool BottomSegmentIterator::hasParseUp() const
{
  return _bottomSegment->hasParseUp();
}

hal_index_t BottomSegmentIterator::getLeftChildIndex(hal_size_t i) const
{
  assert(_startOffset == 0 && _endOffset == 0);
  return _bottomSegment->getLeftChildIndex(i);
}

hal_index_t BottomSegmentIterator::getRightChildIndex(hal_size_t i) const
{
  assert(_startOffset == 0 && _endOffset == 0);
  return _bottomSegment->getRightChildIndex(i);
}

//////////////////////////////////////////////////////////////////////////////
// BOTTOM SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void BottomSegmentIterator::toParent(TopSegmentIteratorPtr ts)
{
  _bottomSegment->setArrayIndex(ts->getGenome()->getParent(),
                                ts->getParentIndex());
  _startOffset = ts->getStartOffset();
  _endOffset = ts->getEndOffset();
  _reversed = ts->getReversed();
  if (ts->getParentReversed() == true)
  {
    toReverse();
  }
  assert (inRange() == true);
}

void 
BottomSegmentIterator::toParseDown(TopSegmentIteratorPtr ts)
{
    Genome* genome = ts->getGenome();
  hal_index_t index = ts->getBottomParseIndex();
  
  _bottomSegment->setArrayIndex(genome, index);
  _reversed = ts->getReversed();
  
  hal_index_t startPos = ts->getStartPosition();
  while (startPos >= _bottomSegment->getStartPosition() + 
         (hal_index_t)_bottomSegment->getLength())
  {
    _bottomSegment->setArrayIndex(genome, ++index);
  }
  
  if (_reversed == false)
  {
    _startOffset = startPos - _bottomSegment->getStartPosition();    
    hal_index_t topEnd = _bottomSegment->getStartPosition() + 
       (hal_index_t)_bottomSegment->getLength();
    hal_index_t botEnd = ts->getStartPosition() + 
       (hal_index_t)ts->getLength();
    _endOffset = max((hal_index_t)0, topEnd - botEnd);
  }
  else
  {
    _startOffset = _bottomSegment->getStartPosition() + 
       _bottomSegment->getLength() - 1
       - startPos;
    hal_index_t topEnd = _bottomSegment->getStartPosition();
    hal_index_t botEnd = ts->getStartPosition() - 
       (hal_index_t)ts->getLength() + 1;
    _endOffset = max((hal_index_t)0, botEnd - topEnd);  
  }
  assert (_startOffset + _endOffset <= _bottomSegment->getLength());
  assert (inRange() == true);
}

BottomSegmentIteratorPtr BottomSegmentIterator::clone() const
{
  assert (inRange() == true);
  BottomSegmentIteratorPtr newIt = 
     getGenome()->getBottomSegmentIterator(getArrayIndex());
  if (_reversed)
  {
    newIt->toReverse();
  }
  newIt->slice(_startOffset, _endOffset);
  return newIt;
}

void BottomSegmentIterator::copy(BottomSegmentIteratorPtr bs)
{
  assert(bs.get() != NULL);
  _bottomSegment->setArrayIndex(bs->getGenome(), bs->getArrayIndex());
  _startOffset = bs->getStartOffset();
  _endOffset = bs->getEndOffset();
  _reversed = bs->getReversed();
}

bool BottomSegmentIterator::equals(BottomSegmentIteratorPtr other) 
  const
{
  assert(_bottomSegment->getGenome() == other->getGenome());
  return getArrayIndex() == other->getArrayIndex();
}

SegmentPtr BottomSegmentIterator::getSegment()
{
  return _bottomSegment;
}

SegmentConstPtr BottomSegmentIterator::getSegment() const
{
  return _bottomSegment;
}





