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
#include "defaultBottomSegmentIterator.h"

using namespace std;
using namespace hal;

DefaultBottomSegmentIterator::
DefaultBottomSegmentIterator(BottomSegment* bottomSegment, 
                             hal_size_t startOffset, 
                             hal_size_t endOffset,
                             bool reversed) :
  DefaultSegmentIterator(startOffset,
                         endOffset,
                         reversed),
  _bottomSegment(bottomSegment)
{

}

DefaultBottomSegmentIterator::~DefaultBottomSegmentIterator()
{

}

SegmentPtr DefaultBottomSegmentIterator::getSegment()
{
  return _bottomSegment;
}

SegmentConstPtr DefaultBottomSegmentIterator::getSegment() const
{
  return _bottomSegment;
}

hal_size_t DefaultBottomSegmentIterator::getNumSegmentsInGenome() const
{
  return getGenome()->getNumBottomSegments();
}
 
//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE OVERRIDE
//////////////////////////////////////////////////////////////////////////////
void DefaultBottomSegmentIterator::print(ostream& os) const
{
  os << "BotSegIt: ";
  DefaultSegmentIterator::print(os);

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
hal_size_t DefaultBottomSegmentIterator::getNumChildren() const
{
  return _bottomSegment->getNumChildren();
}

hal_index_t DefaultBottomSegmentIterator::getChildIndex(hal_size_t i) const
{
  return _bottomSegment->getChildIndex(i);
}

hal_index_t DefaultBottomSegmentIterator::getChildIndexG(
  const Genome* childGenome) const
{
  return _bottomSegment->getChildIndexG(childGenome);
}

bool DefaultBottomSegmentIterator::hasChild(hal_size_t child) const
{
  return _bottomSegment->hasChild(child);
}

bool DefaultBottomSegmentIterator::hasChildG(const Genome* childGenome) const
{
  return _bottomSegment->hasChildG(childGenome);
}

void DefaultBottomSegmentIterator::setChildIndex(hal_size_t i, 
                                              hal_index_t childIndex)
{
  _bottomSegment->setChildIndex(i, childIndex);
}

bool DefaultBottomSegmentIterator::getChildReversed(hal_size_t i) const
{
  return _bottomSegment->getChildReversed(i);
}

void DefaultBottomSegmentIterator::setChildReversed(hal_size_t child,
                                                 bool isReversed)
{
  _bottomSegment->setChildReversed(child, isReversed);
}

hal_index_t DefaultBottomSegmentIterator::getTopParseIndex() const
{
  return _bottomSegment->getTopParseIndex();
}

void DefaultBottomSegmentIterator::setTopParseIndex(hal_index_t parseIndex)
{
  _bottomSegment->setTopParseIndex(parseIndex);
}

hal_offset_t DefaultBottomSegmentIterator::getTopParseOffset() const
{
  return _bottomSegment->getTopParseOffset();
}

bool DefaultBottomSegmentIterator::hasParseUp() const
{
  return _bottomSegment->hasParseUp();
}

hal_index_t DefaultBottomSegmentIterator::getLeftChildIndex(hal_size_t i) const
{
  assert(_startOffset == 0 && _endOffset == 0);
  return _bottomSegment->getLeftChildIndex(i);
}

hal_index_t DefaultBottomSegmentIterator::getRightChildIndex(hal_size_t i) const
{
  assert(_startOffset == 0 && _endOffset == 0);
  return _bottomSegment->getRightChildIndex(i);
}

//////////////////////////////////////////////////////////////////////////////
// BOTTOM SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultBottomSegmentIterator::toParent(TopSegmentIteratorConstPtr ts) const
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
DefaultBottomSegmentIterator::toParseDown(TopSegmentIteratorConstPtr ts) const
{
  const Genome* genome = ts->getGenome();
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

BottomSegmentIteratorPtr DefaultBottomSegmentIterator::copy()
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

BottomSegmentIteratorConstPtr DefaultBottomSegmentIterator::copy() const
{
  assert (inRange() == true);
  BottomSegmentIteratorConstPtr newIt = 
     getGenome()->getBottomSegmentIterator(getArrayIndex());
  if (_reversed)
  {
    newIt->toReverse();
  }
  newIt->slice(_startOffset, _endOffset);
  return newIt;
}

void DefaultBottomSegmentIterator::copy(BottomSegmentIteratorConstPtr bs) const
{
  assert(bs.get() != NULL);
  _bottomSegment->setArrayIndex(bs->getGenome(), bs->getArrayIndex());
  _startOffset = bs->getStartOffset();
  _endOffset = bs->getEndOffset();
  _reversed = bs->getReversed();
}

BottomSegment* DefaultBottomSegmentIterator::getBottomSegment()
{
  // Deprecated now, but so much current code relies on these functions
  return _bottomSegment.get();
}

const BottomSegment* DefaultBottomSegmentIterator::getBottomSegment() const
{
  // Deprecated now, but so much current code relies on these functions
  return _bottomSegment.get();
}

bool DefaultBottomSegmentIterator::equals(BottomSegmentIteratorConstPtr other) 
  const
{
  assert(_bottomSegment->getGenome() == other->getGenome());
  return getArrayIndex() == other->getArrayIndex();
}





