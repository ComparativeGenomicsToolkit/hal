/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hal.h"
#include "defaultSlicedSegment.h"

using namespace std;
using namespace hal;

DefaultSlicedSegment::DefaultSlicedSegment(hal_offset_t startOffset, 
                                           hal_offset_t endOffset,
                                           bool reversed) :
  _startOffset(startOffset),
  _endOffset(endOffset),
  _reversed(reversed)
{

}

DefaultSlicedSegment::~DefaultSlicedSegment()
{

}
   
//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultSlicedSegment::setArrayIndex(Genome* genome, 
                                         hal_index_t arrayIndex)
{
  getSegment()->setArrayIndex(genome, arrayIndex);
}

void DefaultSlicedSegment::setArrayIndex(const Genome* genome, 
                                         hal_index_t arrayIndex) const
{
  getSegment()->setArrayIndex(genome, arrayIndex);
}

const Genome* DefaultSlicedSegment::getGenome() const
{
  return getSegment()->getGenome();
}

Genome* DefaultSlicedSegment::getGenome()
{
  return getSegment()->getGenome();
}

const Sequence* DefaultSlicedSegment::getSequence() const
{
  return getSegment()->getSequence();
}

Sequence* DefaultSlicedSegment::getSequence()
{
  return getSegment()->getSequence();
}

hal_index_t DefaultSlicedSegment::getStartPosition() const
{
  assert (inRange() == true);
  if (_reversed == false)
  {
    return getSegment()->getStartPosition() + _startOffset;
  }
  else
  {
    return getSegment()->getStartPosition() + getSegment()->getLength() - 
       _startOffset - 1;
  }
}

hal_index_t DefaultSlicedSegment::getEndPosition() const
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

hal_size_t DefaultSlicedSegment::getLength() const
{
  assert (inRange() == true);
  return getSegment()->getLength() - _endOffset - _startOffset;
}

void DefaultSlicedSegment::getString(std::string& outString) const
{
  assert (inRange() == true);
  getSegment()->getString(outString);
  if (_reversed == true)
  {
    reverseComplement(outString);
  }
  outString = outString.substr(_startOffset, getLength());
}

void DefaultSlicedSegment::setCoordinates(hal_index_t startPos, 
                                          hal_size_t length)
{
  getSegment()->setCoordinates(startPos, length);
}

hal_index_t DefaultSlicedSegment::getArrayIndex() const
{
  return getSegment()->getArrayIndex();
}

bool DefaultSlicedSegment::leftOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(getSegment()->getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return (hal_index_t)(getStartPosition() + getLength()) <= genomePos;
  }
  else
  {
    return (hal_index_t)getStartPosition() < genomePos;
  }
}

bool DefaultSlicedSegment::rightOf(hal_index_t genomePos) const
{
  assert(genomePos != NULL_INDEX);
  assert(getSegment()->getStartPosition() != NULL_INDEX);
  if (_reversed == false)
  {
    return getStartPosition() > genomePos;
  }
  else
  {
    return getStartPosition() - (hal_index_t)getLength() >= genomePos;
  }
}

bool DefaultSlicedSegment::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

bool DefaultSlicedSegment::isFirst() const
{
  return !_reversed ? getSegment()->isFirst() : getSegment()->isLast();
}

bool DefaultSlicedSegment::isLast() const
{
  return !_reversed ? getSegment()->isLast() : getSegment()->isFirst();
}

bool DefaultSlicedSegment::isMissingData(double nThreshold) const
{
  return getSegment()->isMissingData(nThreshold);
}

bool DefaultSlicedSegment::isTop() const
{
  return getSegment()->isTop();
}

hal_size_t DefaultSlicedSegment::getMappedSegments(
  const Genome* tgtGenome,
  std::vector<MappedSegmentConstPtr>& outSegments,
  bool doDupes) const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////
// SLICED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultSlicedSegment::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
//  swap(_startOffset, _endOffset);
}

hal_offset_t DefaultSlicedSegment::getStartOffset() const
{
  return _startOffset;
}

hal_offset_t DefaultSlicedSegment::getEndOffset() const
{
  return _endOffset;
}

void DefaultSlicedSegment::slice(hal_offset_t startOffset, 
                                   hal_offset_t endOffset) const
{
  assert(startOffset < getSegment()->getLength());
  assert(endOffset < getSegment()->getLength());
  _startOffset = startOffset;
  _endOffset = endOffset;
}


bool DefaultSlicedSegment::getReversed() const
{
  return _reversed;
}

