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

void HDF5TopSegmentIterator::toLeft() const
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
}

void HDF5TopSegmentIterator::toRight() const
{
  if (_endOffset == 0)
  {
    ++_topSegment._index;
    _startOffset = 0;
  }
  else
  {
    _startOffset =  _topSegment.getLength() - _endOffset - 1;
    _endOffset = 0;
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
  _startOffset = h5bs->_bottomSegment.getTopParseOffset();
  _endOffset = min(static_cast<hal_size_t>(h5bs->_endOffset),
                   _topSegment.getLength());
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
  assert (inRange() == true);
  return &_topSegment;
}

const TopSegment* HDF5TopSegmentIterator::getTopSegment() const
{
  assert (inRange() == true);
  return &_topSegment;
}
