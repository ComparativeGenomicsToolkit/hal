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

void HDF5BottomSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  if (_reversed == false)
  {
    if (_startOffset == 0)
    {
      --_bottomSegment._index;
      _endOffset = 0;
    }
    else
    {
      _endOffset = _bottomSegment.getLength() - _startOffset;
      _startOffset = 0;
    }
    if (leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      assert(_bottomSegment.getStartPosition() <= leftCutoff);
      _startOffset = leftCutoff - _bottomSegment.getStartPosition();
    }
  }
  
  else
  {
    if (_startOffset == 0)
    {
      ++_bottomSegment._index;
      _endOffset = 0;
    }
    else
    {
      _endOffset = _bottomSegment.getLength() - _startOffset;
      _startOffset = 0;
    }    
    if ((hal_size_t)_bottomSegment._index < _bottomSegment._array->getSize() &&
        leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      _startOffset = _bottomSegment.getStartPosition() + 
         _bottomSegment.getLength() - 1 - leftCutoff;
    }
  }
   assert((hal_size_t)_bottomSegment._index >= 
          _bottomSegment._array->getSize() ||
          _bottomSegment._index < 0 || 
          _startOffset + _endOffset <= _bottomSegment.getLength());
}

void HDF5BottomSegmentIterator::toRight(hal_index_t rightCutoff) const  
{
  if (_reversed == false)
  {
    if (_endOffset == 0)
    {
      ++_bottomSegment._index;
      _startOffset = 0;
    }
    else
    {
      _startOffset =  _bottomSegment.getLength() - _endOffset;
      _endOffset = 0;
    }
    
    if ((hal_size_t)_bottomSegment._index < _bottomSegment._array->getSize() &&
        rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = _bottomSegment.getStartPosition() +
         _bottomSegment.getLength() - rightCutoff - 1;
    }
  }
  
  else
  {
    if (_endOffset == 0)
    {
      --_bottomSegment._index;
      _startOffset = 0;
    }
    else
    {
      _startOffset =  _bottomSegment.getLength() - _endOffset;
      _endOffset = 0;
    }

    if (rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = rightCutoff - _bottomSegment.getStartPosition(); 
    }
  }
  assert ((hal_size_t)_bottomSegment._index >= 
          _bottomSegment._array->getSize() ||
          _bottomSegment._index < 0 || 
          _startOffset + _endOffset <= _bottomSegment.getLength());
}

void HDF5BottomSegmentIterator::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
  swap(_startOffset, _endOffset);
}

void HDF5BottomSegmentIterator::toSite(hal_index_t position, bool slice) const
{
  hal_index_t len = (hal_index_t)_bottomSegment._genome->getSequenceLength();
  hal_index_t nseg = (hal_index_t)_bottomSegment._genome->getNumBottomSegments();

  assert(len != 0);
  double avgLen = (double)len / (double)nseg;
  hal_index_t hint = (hal_index_t)
     min(nseg - 1., avgLen * ((double)position / (double)len));
  _bottomSegment._index = hint;
  _startOffset = 0;
  _endOffset = 0;

  // out of range
  if (position < 0)
  {
    _bottomSegment._index = -1;
    return;
  }
  else if (position >= len)
  {
    _bottomSegment._index = len;
    return;
  }
  
  hal_index_t left = 0;
  hal_index_t right = nseg - 1;
  assert(_bottomSegment._index  >= 0 &&  _bottomSegment._index < nseg);
  
  while (overlaps(position) == false)
  {
    assert(left != right);
    if (rightOf(position) == true)
    {
      right = _bottomSegment._index;
      _bottomSegment._index -= max((_bottomSegment._index - left) / 2,
                                   (hal_index_t)1);
      assert(_bottomSegment._index  >= 0 &&  _bottomSegment._index < nseg);
    }
    else
    {
      assert(leftOf(position) == true);
      left = _bottomSegment._index;
      _bottomSegment._index += max((right - _bottomSegment._index) / 2,
                                   (hal_index_t)1);
      assert(_bottomSegment._index  >= 0 &&  _bottomSegment._index < nseg);
    }
  }
  
  assert(overlaps(position) == true);
  
  if (slice == true)
  {
     _startOffset = position - _bottomSegment.getStartPosition();
     _endOffset = _bottomSegment.getStartPosition() + _bottomSegment.getLength()
        - position - 1;
     if (_reversed)
     {
       swap(_startOffset, _endOffset);
     }
  }  
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

void HDF5BottomSegmentIterator::slice(hal_offset_t startOffset, 
                                      hal_offset_t endOffset) const
{
  assert(startOffset < _bottomSegment.getLength());
  assert(endOffset < _bottomSegment.getLength());
  _startOffset = startOffset;
  _endOffset = endOffset;
}

hal_index_t HDF5BottomSegmentIterator::getStartPosition() const
{
  assert (inRange() == true);
  if (_reversed == false)
  {
    return _bottomSegment.getStartPosition() + _startOffset;
  }
  else
  {
    return _bottomSegment.getStartPosition() + _bottomSegment.getLength() - 
       _startOffset - 1;
  }
}


hal_size_t HDF5BottomSegmentIterator::getLength() const
{
  assert (inRange() == true);
  assert (_endOffset + _startOffset <= _bottomSegment.getLength());
  return _bottomSegment.getLength() - _endOffset - _startOffset;
}

hal_bool_t HDF5BottomSegmentIterator::getReversed() const
{
  return _reversed;
}

void HDF5BottomSegmentIterator::getString(std::string& outString) const
{
  assert (inRange() == true);
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_bottomSegment._genome), 
                     getStartPosition());
  if (_reversed == true)
  {
    di.toReverse();
  }
  di.readString(outString, getLength()); 
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
  if (h5ts->_topSegment.getParentReversed() == true)
  {
    toReverse();
  }
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
  _reversed = h5ts->_reversed;
  
  hal_index_t startPos = h5ts->getStartPosition();
  while (startPos >= _bottomSegment.getStartPosition() + 
         (hal_index_t)_bottomSegment.getLength())
  {
    ++_bottomSegment._index;
  }
  
  if (_reversed == false)
  {
    _startOffset = startPos - _bottomSegment.getStartPosition();    
    hal_index_t topEnd = _bottomSegment.getStartPosition() + 
       (hal_index_t)_bottomSegment.getLength();
    hal_index_t botEnd = h5ts->getStartPosition() + 
       (hal_index_t)h5ts->getLength();
    _endOffset = max((hal_index_t)0, topEnd - botEnd);
  }
  else
  {
    _startOffset = _bottomSegment.getStartPosition() + 
       _bottomSegment.getLength() - 1
       - startPos;
    hal_index_t topEnd = _bottomSegment.getStartPosition();
    hal_index_t botEnd = h5ts->getStartPosition() - 
       (hal_index_t)h5ts->getLength() + 1;
    _endOffset = max((hal_index_t)0, botEnd - topEnd);  
  }
  assert (_startOffset + _endOffset <= _bottomSegment.getLength());
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

void HDF5BottomSegmentIterator::copy(BottomSegmentIteratorConstPtr bs) const
{
  assert(bs.get() != NULL);
  const HDF5BottomSegmentIterator* h5bs = 
     reinterpret_cast<const HDF5BottomSegmentIterator*>(bs.get());
  _bottomSegment._array = h5bs->_bottomSegment._array;
  _bottomSegment._genome = h5bs->_bottomSegment._genome;
  _bottomSegment._index = h5bs->_bottomSegment._index;
  _startOffset = h5bs->_startOffset;
  _endOffset = h5bs->_endOffset;
  _reversed = h5bs->_reversed;
}

BottomSegment* HDF5BottomSegmentIterator::getBottomSegment()
{
  return &_bottomSegment;
}

const BottomSegment* HDF5BottomSegmentIterator::getBottomSegment() const
{
  return &_bottomSegment;
}

bool HDF5BottomSegmentIterator::equals(BottomSegmentIteratorConstPtr other) const
{
  const HDF5BottomSegmentIterator* h5Other = reinterpret_cast<
     const HDF5BottomSegmentIterator*>(other.get());
  assert(_bottomSegment.getGenome() == h5Other->_bottomSegment.getGenome());
  return _bottomSegment._index == h5Other->_bottomSegment._index;
}

bool HDF5BottomSegmentIterator::hasChild(hal_size_t child) const
{
  assert (inRange() == true);
  return _bottomSegment.getChildIndex(child) != NULL_INDEX;
}

bool HDF5BottomSegmentIterator::hasParseUp() const
{
  assert (inRange() == true);
  assert (_bottomSegment.getTopParseIndex() == NULL_INDEX || 
          _bottomSegment.getGenome()->getParent() != NULL);
  return _bottomSegment.getTopParseIndex() != NULL_INDEX;
}

bool HDF5BottomSegmentIterator::leftOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_bottomSegment.getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return (hal_index_t)(getStartPosition() + getLength()) <= genomePos;
  }
  else
  {
    return (hal_index_t)getStartPosition() < genomePos;
  }
}

bool HDF5BottomSegmentIterator::rightOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_bottomSegment.getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return getStartPosition() > genomePos;
  }
  else
  {
    return getStartPosition() - (hal_index_t)getLength() >= genomePos;
  }
}

bool HDF5BottomSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

bool HDF5BottomSegmentIterator::isFirst() const
{
  return !_reversed ? _bottomSegment.isFirst() : _bottomSegment.isLast();
}

bool HDF5BottomSegmentIterator::isLast() const
{
  return !_reversed ? _bottomSegment.isLast() : _bottomSegment.isFirst();
}



