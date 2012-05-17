/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hdf5TopSegmentIterator.h"
#include "hdf5BottomSegmentIterator.h"
#include "hdf5DNAIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

HDF5TopSegmentIterator::HDF5TopSegmentIterator(HDF5Genome* genome, 
                                               hal_index_t index,
                                               hal_offset_t startOffset, 
                                               hal_offset_t endOffset,
                                               hal_bool_t reversed) :
  _topSegment(genome, &genome->_topArray, index),
  _startOffset(startOffset),
  _endOffset(endOffset),
  _reversed(reversed)
{

}

HDF5TopSegmentIterator::~HDF5TopSegmentIterator()
{
}
   
// ITERATOR METHODS

void HDF5TopSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  if (_startOffset == 0)
  {
    --_topSegment._index;
    _endOffset = 0;
  }
  else
  {
    _endOffset = _topSegment.getLength() - _startOffset - 1;
    _startOffset = 0;
  }
  if (leftCutoff != NULL_INDEX && overlaps(leftCutoff))
  {
    _startOffset = leftCutoff - _topSegment.getStartPosition();
  }
}

void HDF5TopSegmentIterator::toRight(hal_index_t rightCutoff) const
{
  if (_endOffset == 0)
  {
    ++_topSegment._index;
    _startOffset = 0;
  }
  else
  {
    _startOffset =  _topSegment.getLength() - _endOffset;
    _endOffset = 0;
  }

  if ((hal_size_t)_topSegment._index < _topSegment._array->getSize() &&
      rightCutoff != NULL_INDEX && overlaps(rightCutoff))
  {
    _endOffset = _topSegment.getStartPosition() + 
       _topSegment.getLength() - rightCutoff - 1;
  }
}

void HDF5TopSegmentIterator::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
  swap(_startOffset, _endOffset);
}

void HDF5TopSegmentIterator::toNextParalogy() const
{
  assert(_topSegment.getNextParalogyIndex() != NULL_INDEX);
  _topSegment._index = _topSegment.getNextParalogyIndex();
/*  if(_topSegment.getGetParalogyReversed() == true)
    {
    _topSegment._revrsed = !_topSegment._reversed;
    }
*/
}

hal_offset_t HDF5TopSegmentIterator::getStartOffset() const
{
  assert (inRange() == true);
  return _startOffset;
}

hal_offset_t HDF5TopSegmentIterator::getEndOffset() const
{
  assert (inRange() == true);
  return _endOffset;
}

void HDF5TopSegmentIterator::slice(hal_offset_t startOffset, 
                                   hal_offset_t endOffset) const
{
  assert(startOffset < _topSegment.getLength());
  assert(endOffset < _topSegment.getLength());
  _startOffset = startOffset;
  _endOffset = endOffset;
}

hal_index_t HDF5TopSegmentIterator::getStartPosition() const
{
  assert (inRange() == true);
  return _topSegment.getStartPosition() + _startOffset;
}

hal_size_t HDF5TopSegmentIterator::getLength() const
{
  assert (inRange() == true);
  return _topSegment.getLength() - _endOffset - _startOffset;
}

hal_bool_t HDF5TopSegmentIterator::getReversed() const
{
  assert (inRange() == true);
  return _reversed;
}

void HDF5TopSegmentIterator::getString(std::string& outString) const
{
  assert (inRange() == true);
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_topSegment._genome), 
                     _topSegment.getStartPosition() + _startOffset);
  di.readString(outString, getLength(), _reversed); 
}

//TOP ITERATOR METHODS


void HDF5TopSegmentIterator::toChild(BottomSegmentIteratorConstPtr bs, 
                                     hal_size_t child) const
{
  assert (inRange() == true);
  const HDF5BottomSegmentIterator* h5bs = 
     reinterpret_cast<const HDF5BottomSegmentIterator*>(bs.get());
  _topSegment._genome = static_cast<HDF5Genome*>(
    h5bs->_bottomSegment._genome->getChild(child));
  _topSegment._index = h5bs->_bottomSegment.getChildIndex(child);
  _topSegment._array = &_topSegment._genome->_topArray;
  _startOffset = h5bs->_startOffset;
  _endOffset = h5bs->_endOffset;
  _reversed = h5bs->_reversed;
}

void HDF5TopSegmentIterator::toParseUp(BottomSegmentIteratorConstPtr bs) const
{
  assert (inRange() == true);
  const HDF5BottomSegmentIterator* h5bs = 
     reinterpret_cast<const HDF5BottomSegmentIterator*>(bs.get());
  
  _topSegment._genome = h5bs->_bottomSegment._genome;
  _topSegment._index = h5bs->_bottomSegment.getTopParseIndex();
  _topSegment._array = &h5bs->_bottomSegment._genome->_topArray;
  
  hal_index_t startPos = h5bs->getStartPosition();
  hal_index_t endPos = startPos + h5bs->getLength() - 1;

  while (startPos >= _topSegment.getStartPosition() + 
         (hal_index_t)_topSegment.getLength())
  {
    ++_topSegment._index;
  }
  
  _startOffset = startPos - _topSegment.getStartPosition();

  hal_index_t newEndPos = _topSegment.getStartPosition() +
     _topSegment.getLength() - 1;
  _endOffset = max((hal_index_t)0, newEndPos - endPos);

  _reversed = h5bs->_reversed;
}

TopSegmentIteratorPtr HDF5TopSegmentIterator::copy()
{
  HDF5TopSegmentIterator* newIt = 
     new HDF5TopSegmentIterator(_topSegment._genome, 
                                _topSegment._index,
                                _startOffset,
                                _endOffset,
                                _reversed);
  return TopSegmentIteratorPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5TopSegmentIterator::copy() const
{
  const HDF5TopSegmentIterator* newIt = 
     new HDF5TopSegmentIterator(_topSegment._genome, 
                                _topSegment._index,
                                _startOffset,
                                _endOffset,
                                _reversed);
  return TopSegmentIteratorConstPtr(newIt);
}

TopSegment* HDF5TopSegmentIterator::getTopSegment()
{
  return &_topSegment;
}

const TopSegment* HDF5TopSegmentIterator::getTopSegment() const
{
  return &_topSegment;
}

bool HDF5TopSegmentIterator::equals(TopSegmentIteratorConstPtr other) const
{
  const HDF5TopSegmentIterator* h5Other = reinterpret_cast<
     const HDF5TopSegmentIterator*>(other.get());
  assert(_topSegment.getGenome() == h5Other->_topSegment.getGenome());
  return _topSegment._index == h5Other->_topSegment._index;
}

bool HDF5TopSegmentIterator::hasParent() const
{
  assert(inRange() == true);
  return _topSegment.getParentIndex() != NULL_INDEX;
}

bool HDF5TopSegmentIterator::hasParseDown() const
{
  assert (inRange() == true);
  assert (_topSegment.getBottomParseIndex() == NULL_INDEX ||
          _topSegment.getGenome()->getNumChildren() > 0);
  return _topSegment.getBottomParseIndex() != NULL_INDEX;
}

bool HDF5TopSegmentIterator::leftOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_topSegment.getStartPosition() != NULL_INDEX);
  return (hal_index_t)(_topSegment.getStartPosition() + _startOffset + 
                       getLength()) < genomePos;
}

bool HDF5TopSegmentIterator::rightOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_topSegment.getStartPosition() != NULL_INDEX);
  return (hal_index_t)(_topSegment.getStartPosition() + _startOffset) >
     genomePos;
}

bool HDF5TopSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}
