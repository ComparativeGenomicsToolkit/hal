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
  _bottomSegment(bottomSegment),
  _startOffset(startOffset),
  _endOffset(endOffset),
  _reversed(reversed)
{

}

DefaultBottomSegmentIterator::~DefaultBottomSegmentIterator()
{
}
   
//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultBottomSegmentIterator::setArrayIndex( Genome* genome, 
                                                  hal_index_t arrayIndex)
{
  _bottomSegment->setArrayIndex(genome, arrayIndex);
}

void DefaultBottomSegmentIterator::setArrayIndex(const Genome* genome, 
                                                 hal_index_t arrayIndex) const
{
  _bottomSegment->setArrayIndex(genome, arrayIndex);
}

const Genome* DefaultBottomSegmentIterator::getGenome() const
{
  return _bottomSegment->getGenome();
}

Genome* DefaultBottomSegmentIterator::getGenome()
{
  return _bottomSegment->getGenome();
}

const Sequence* DefaultBottomSegmentIterator::getSequence() const
{
  return _bottomSegment->getSequence();
}

Sequence* DefaultBottomSegmentIterator::getSequence()
{
  return _bottomSegment->getSequence();
}

hal_index_t DefaultBottomSegmentIterator::getStartPosition() const
{
  assert (inRange() == true);
  if (_reversed == false)
  {
    return _bottomSegment->getStartPosition() + _startOffset;
  }
  else
  {
    return _bottomSegment->getStartPosition() + _bottomSegment->getLength() - 
       _startOffset - 1;
  }
}

hal_index_t DefaultBottomSegmentIterator::getEndPosition() const
{
  assert (inRange() == true);
  if (_reversed == false)
  {
    return getStartPosition() + (hal_index_t)(getLength() - 1);
  }
  else
  {
    return getStartPosition() - (hal_index_t)(getLength() - 1);
  }
}

hal_size_t DefaultBottomSegmentIterator::getLength() const
{
  assert (inRange() == true);
  assert (_endOffset + _startOffset <= _bottomSegment->getLength());
  return _bottomSegment->getLength() - _endOffset - _startOffset;
}

void DefaultBottomSegmentIterator::getString(std::string& outString) const
{
  assert (inRange() == true);
  _bottomSegment->getString(outString);
  if (_reversed == true)
  {
    reverseComplement(outString);
  }
  outString = outString.substr(_startOffset, getLength());
}

void DefaultBottomSegmentIterator::setCoordinates(hal_index_t startPos, 
                                               hal_size_t length)
{
  _bottomSegment->setCoordinates(startPos, length);
}

hal_index_t DefaultBottomSegmentIterator::getArrayIndex() const
{
  return _bottomSegment->getArrayIndex();
}

bool DefaultBottomSegmentIterator::leftOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_bottomSegment->getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return (hal_index_t)(getStartPosition() + getLength()) <= genomePos;
  }
  else
  {
    return (hal_index_t)getStartPosition() < genomePos;
  }
}

bool DefaultBottomSegmentIterator::rightOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(_bottomSegment->getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return getStartPosition() > genomePos;
  }
  else
  {
    return getStartPosition() - (hal_index_t)getLength() >= genomePos;
  }
}

bool DefaultBottomSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

bool DefaultBottomSegmentIterator::isFirst() const
{
  return !_reversed ? _bottomSegment->isFirst() : _bottomSegment->isLast();
}

bool DefaultBottomSegmentIterator::isLast() const
{
  return !_reversed ? _bottomSegment->isLast() : _bottomSegment->isFirst();
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
// SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultBottomSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  if (_reversed == false)
  {
    if (_startOffset == 0)
    {
      _bottomSegment->setArrayIndex(getGenome(), getArrayIndex() - 1);
      _endOffset = 0;
    }
    else
    {
      _endOffset = _bottomSegment->getLength() - _startOffset;
      _startOffset = 0;
    }
    if (leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      assert(_bottomSegment->getStartPosition() <= leftCutoff);
      _startOffset = leftCutoff - _bottomSegment->getStartPosition();
    }
  }
  
  else
  {
    if (_startOffset == 0)
    {
      _bottomSegment->setArrayIndex(getGenome(), getArrayIndex() + 1);
      _endOffset = 0;
    }
    else
    {
      _endOffset = _bottomSegment->getLength() - _startOffset;
      _startOffset = 0;
    }    
    if ((hal_size_t)_bottomSegment->getArrayIndex() <
        _bottomSegment->getGenome()->getNumBottomSegments() &&
        leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      _startOffset = _bottomSegment->getStartPosition() + 
         _bottomSegment->getLength() - 1 - leftCutoff;
    }
  }
  assert((hal_size_t)_bottomSegment->getArrayIndex() >= 
         _bottomSegment->getGenome()->getNumBottomSegments() ||
         _bottomSegment->getArrayIndex() < 0 || 
          _startOffset + _endOffset <= _bottomSegment->getLength());
}

void DefaultBottomSegmentIterator::toRight(hal_index_t rightCutoff) const  
{
  if (_reversed == false)
  {
    if (_endOffset == 0)
    {
      _bottomSegment->setArrayIndex(getGenome(), getArrayIndex() + 1);
      _startOffset = 0;
    }
    else
    {
      _startOffset =  _bottomSegment->getLength() - _endOffset;
      _endOffset = 0;
    }
    
    if ((hal_size_t)_bottomSegment->getArrayIndex() < 
        _bottomSegment->getGenome()->getNumBottomSegments() &&
        rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = _bottomSegment->getStartPosition() +
         _bottomSegment->getLength() - rightCutoff - 1;
    }
  }
  
  else
  {
    if (_endOffset == 0)
    {
      _bottomSegment->setArrayIndex(getGenome(), getArrayIndex() - 1);
      _startOffset = 0;
    }
    else
    {
      _startOffset =  _bottomSegment->getLength() - _endOffset;
      _endOffset = 0;
    }

    if (rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = rightCutoff - _bottomSegment->getStartPosition(); 
    }
  }
  assert ((hal_size_t)_bottomSegment->getArrayIndex() >= 
          _bottomSegment->getGenome()->getNumBottomSegments() ||
          _bottomSegment->getArrayIndex() < 0 || 
          _startOffset + _endOffset <= _bottomSegment->getLength());
}

void DefaultBottomSegmentIterator::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
  swap(_startOffset, _endOffset);
}

void DefaultBottomSegmentIterator::toSite(hal_index_t position, bool slice) const
{
  const Genome* genome = getGenome();
  hal_index_t len = (hal_index_t)genome->getSequenceLength();
  hal_index_t nseg = (hal_index_t)genome->getNumBottomSegments();

  assert(len != 0);
  double avgLen = (double)len / (double)nseg;
  hal_index_t hint = (hal_index_t)
     min(nseg - 1., avgLen * ((double)position / (double)len));
  _bottomSegment->setArrayIndex(genome, hint);
  _startOffset = 0;
  _endOffset = 0;

  // out of range
  if (position < 0)
  {
    _bottomSegment->setArrayIndex(genome, NULL_INDEX);
    return;
  }
  else if (position >= len)
  {
    _bottomSegment->setArrayIndex(genome, len);
    return;
  }
  
  hal_index_t left = 0;
  hal_index_t right = nseg - 1;
  assert(_bottomSegment->getArrayIndex() >= 0 &&  
         _bottomSegment->getArrayIndex() < nseg);
  
  while (overlaps(position) == false)
  {
    assert(left != right);
    if (rightOf(position) == true)
    {
      right = _bottomSegment->getArrayIndex();
      hal_index_t delta =  max((_bottomSegment->getArrayIndex() - left) / 2,
                               (hal_index_t)1);
      _bottomSegment->setArrayIndex(genome, 
                                    _bottomSegment->getArrayIndex() - delta);
      assert(_bottomSegment->getArrayIndex()  >= 0 &&  
             _bottomSegment->getArrayIndex() < nseg);
    }
    else
    {
      assert(leftOf(position) == true);
      left = _bottomSegment->getArrayIndex();
      hal_index_t delta = max((right - _bottomSegment->getArrayIndex()) / 2,
                              (hal_index_t)1);
      _bottomSegment->setArrayIndex(genome, 
                                    _bottomSegment->getArrayIndex() + delta);
      assert(_bottomSegment->getArrayIndex()  >= 0 &&
             _bottomSegment->getArrayIndex() < nseg);
    }
  }
  
  assert(overlaps(position) == true);
  
  if (slice == true)
  {
     _startOffset = position - _bottomSegment->getStartPosition();
     _endOffset = _bottomSegment->getStartPosition() + _bottomSegment->getLength()
        - position - 1;
     if (_reversed)
     {
       swap(_startOffset, _endOffset);
     }
  }  
}

hal_offset_t DefaultBottomSegmentIterator::getStartOffset() const
{
  return _startOffset;
}

hal_offset_t DefaultBottomSegmentIterator::getEndOffset() const
{
  return _endOffset;
}

void DefaultBottomSegmentIterator::slice(hal_offset_t startOffset, 
                                      hal_offset_t endOffset) const
{
  assert(startOffset < _bottomSegment->getLength());
  assert(endOffset < _bottomSegment->getLength());
  _startOffset = startOffset;
  _endOffset = endOffset;
}

bool DefaultBottomSegmentIterator::getReversed() const
{
  return _reversed;
}


//////////////////////////////////////////////////////////////////////////////
// BOTTOM SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultBottomSegmentIterator::toParent(TopSegmentIteratorConstPtr ts) const
{
  assert (inRange() == true);
  _bottomSegment->setArrayIndex(ts->getGenome()->getParent(),
                                ts->getParentIndex());
  _startOffset = ts->getStartOffset();
  _endOffset = ts->getEndOffset();
  _reversed = ts->getReversed();
  if (ts->getParentReversed() == true)
  {
    toReverse();
  }
}

void 
DefaultBottomSegmentIterator::toParseDown(TopSegmentIteratorConstPtr ts) const
{
  assert (inRange() == true);

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





