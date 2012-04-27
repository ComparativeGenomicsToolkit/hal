/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hdf5BottomSegmentIterator.h"
#include "hdf5TopSegmentIterator.h"
#include "hdf5DNAIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

HDF5BottomSegmentIterator::HDF5BottomSegmentIterator(HDF5Genome* genome, 
                                               hal_index_t index,
                                               hal_size_t startOffset, 
                                               hal_size_t endOffset,
                                               hal_bool_t reversed) :
  _bottomSegment(genome, &genome->_bottomArray, index),
  _startOffset(startOffset),
  _endOffset(endOffset),
  _reversed(reversed)
{
 
}

HDF5BottomSegmentIterator::~HDF5BottomSegmentIterator()
{
}
   
// ITERATOR METHODS

void HDF5BottomSegmentIterator::toLeft() const
{
   if (_startOffset == 0)
  {
    --_bottomSegment._index;
    _endOffset = 0;
  }
  else
  {
    _endOffset = _bottomSegment.getLength() - _startOffset - 1;
    _startOffset = 0;
  }
}

void HDF5BottomSegmentIterator::toRight() const
{
   if (_endOffset == 0)
  {
    ++_bottomSegment._index;
    _startOffset = 0;
  }
  else
  {
    _startOffset =  _bottomSegment.getLength() - _endOffset - 1;
    _endOffset = 0;
  }
}

void HDF5BottomSegmentIterator::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
  swap(_startOffset, _endOffset);
}

void HDF5BottomSegmentIterator::toNextParalogy() const
{
  assert (inRange() == true);
  assert(_bottomSegment.getNextParalogyIndex() != NULL_INDEX);
  _bottomSegment._index = _bottomSegment.getNextParalogyIndex();
/*  if(_bottomSegment.getGetParalogyReversed() == true)
  {
  _bottomSegment._revrsed = !_bottomSegment._reversed;
}
*/
}

hal_offset_t HDF5BottomSegmentIterator::getStartOffset() const
{
  assert (inRange() == true);
  return _startOffset;
}

hal_offset_t HDF5BottomSegmentIterator::getEndOffset() const
{
  assert (inRange() == true);
  return _endOffset;
}

hal_size_t HDF5BottomSegmentIterator::getLength() const
{
  assert (inRange() == true);
  return _bottomSegment.getLength() - _endOffset - _startOffset;
}

hal_bool_t HDF5BottomSegmentIterator::getReversed() const
{
  assert (inRange() == true);
  return _reversed;
}

void HDF5BottomSegmentIterator::getString(std::string& outString) const
{
  assert (inRange() == true);
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_bottomSegment._genome), 
                     _bottomSegment.getStartPosition() + _startOffset);
  di.readString(outString, getLength(), _reversed); 
}

//BOTTOM ITERATOR METHODS


void HDF5BottomSegmentIterator::toParent(TopSegmentIteratorConstPtr ts) const
{
  assert (inRange() == true);
  const HDF5TopSegmentIterator* h5ts = 
     reinterpret_cast<const HDF5TopSegmentIterator*>(ts.get());
  _bottomSegment._genome = static_cast<HDF5Genome*>(
    h5ts->_topSegment._genome->getParent());
  _bottomSegment._index = h5ts->_topSegment.getParentIndex();
  _bottomSegment._array = &_bottomSegment._genome->_bottomArray;
  _startOffset = h5ts->_startOffset;
  _endOffset = h5ts->_endOffset;
  _reversed = h5ts->_reversed;
}

void 
HDF5BottomSegmentIterator::toParseDown(TopSegmentIteratorConstPtr ts) const
{
  assert (inRange() == true);
  const HDF5TopSegmentIterator* h5ts = 
     reinterpret_cast<const HDF5TopSegmentIterator*>(ts.get());
  _bottomSegment._genome = h5ts->_topSegment._genome;
  _bottomSegment._index = h5ts->_topSegment.getBottomParseIndex();
  _bottomSegment._array = &h5ts->_topSegment._genome->_bottomArray;
  _startOffset = h5ts->_topSegment.getBottomParseOffset();
  _endOffset = min(static_cast<hal_size_t>(h5ts->_endOffset),
                   _bottomSegment.getLength());
  _reversed = h5ts->_reversed;
}

BottomSegmentIteratorPtr HDF5BottomSegmentIterator::copy()
{
  assert (inRange() == true);
  HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(_bottomSegment._genome, 
                                _bottomSegment._index,
                                _startOffset,
                                _endOffset,
                                _reversed);
  return BottomSegmentIteratorPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5BottomSegmentIterator::copy() const
{
  assert (inRange() == true);
  const HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(_bottomSegment._genome, 
                                _bottomSegment._index,
                                _startOffset,
                                _endOffset,
                                _reversed);
  return BottomSegmentIteratorConstPtr(newIt);
}

BottomSegment* HDF5BottomSegmentIterator::getBottomSegment()
{
  assert (inRange() == true);
  return &_bottomSegment;
}

const BottomSegment* HDF5BottomSegmentIterator::getBottomSegment() const
{
  assert (inRange() == true);
  return &_bottomSegment;
}

bool HDF5BottomSegmentIterator::equals(BottomSegmentIteratorConstPtr other) const
{
  const HDF5BottomSegmentIterator* h5Other = reinterpret_cast<
     const HDF5BottomSegmentIterator*>(other.get());
  assert(_bottomSegment.getGenome() == h5Other->_bottomSegment.getGenome());
  return _bottomSegment._index == h5Other->_bottomSegment._index;
}
