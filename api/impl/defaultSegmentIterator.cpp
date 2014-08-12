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
#include "defaultSegmentIterator.h"
#include "defaultMappedSegment.h"

using namespace std;
using namespace hal;

DefaultSegmentIterator::DefaultSegmentIterator(hal_offset_t startOffset, 
                                               hal_offset_t endOffset,
                                               bool reversed) :
  _startOffset(startOffset),
  _endOffset(endOffset),
  _reversed(reversed)
{

}

DefaultSegmentIterator::~DefaultSegmentIterator()
{

}

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultSegmentIterator::setArrayIndex(Genome* genome, 
                                         hal_index_t arrayIndex)
{
  getSegment()->setArrayIndex(genome, arrayIndex);
}

void DefaultSegmentIterator::setArrayIndex(const Genome* genome, 
                                         hal_index_t arrayIndex) const
{
  getSegment()->setArrayIndex(genome, arrayIndex);
}

const Genome* DefaultSegmentIterator::getGenome() const
{
  return getSegment()->getGenome();
}

Genome* DefaultSegmentIterator::getGenome()
{
  return getSegment()->getGenome();
}

const Sequence* DefaultSegmentIterator::getSequence() const
{
  return getSegment()->getSequence();
}

Sequence* DefaultSegmentIterator::getSequence()
{
  return getSegment()->getSequence();
}

hal_index_t DefaultSegmentIterator::getStartPosition() const
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

hal_index_t DefaultSegmentIterator::getEndPosition() const
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

hal_size_t DefaultSegmentIterator::getLength() const
{
  assert (inRange() == true);
  return getSegment()->getLength() - _endOffset - _startOffset;
}

void DefaultSegmentIterator::getString(std::string& outString) const
{
  assert (inRange() == true);
  getSegment()->getString(outString);
  if (_reversed == true)
  {
    reverseComplement(outString);
  }
  outString = outString.substr(_startOffset, getLength());
}

void DefaultSegmentIterator::setCoordinates(hal_index_t startPos, 
                                          hal_size_t length)
{
  getSegment()->setCoordinates(startPos, length);
}

hal_index_t DefaultSegmentIterator::getArrayIndex() const
{
  return getSegment()->getArrayIndex();
}

bool DefaultSegmentIterator::leftOf(hal_index_t genomePos) const
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

bool DefaultSegmentIterator::rightOf(hal_index_t genomePos) const
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

bool DefaultSegmentIterator::overlaps(hal_index_t genomePos) const
{
  return !leftOf(genomePos) && !rightOf(genomePos);
}

bool DefaultSegmentIterator::isFirst() const
{
  return !_reversed ? getSegment()->isFirst() : getSegment()->isLast();
}

bool DefaultSegmentIterator::isLast() const
{
  return !_reversed ? getSegment()->isLast() : getSegment()->isFirst();
}

bool DefaultSegmentIterator::isMissingData(double nThreshold) const
{
  return getSegment()->isMissingData(nThreshold);
}

bool DefaultSegmentIterator::isTop() const
{
  return getSegment()->isTop();
}

hal_size_t DefaultSegmentIterator::getMappedSegments(
  set<MappedSegmentConstPtr>& outSegments,
  const Genome* tgtGenome,
  const set<const Genome*>* genomesOnPath,
  bool doDupes,
  hal_size_t minLength,
  const Genome *coalescenceLimit,
  const Genome *mrca) const
{
  assert(tgtGenome != NULL);

  if (mrca == NULL)
  {
    set<const Genome*> inputSet;
    inputSet.insert(getGenome());
    inputSet.insert(tgtGenome);
    mrca = getLowestCommonAncestor(inputSet);
  }

  if (coalescenceLimit == NULL)
  {
    coalescenceLimit = mrca;
  }

  // Get the path from the coalescence limit to the target (necessary
  // for choosing which children to move through to get to the
  // target).
  set<const Genome*> pathSet;
  if (genomesOnPath == NULL)
  {
    set<const Genome*> inputSet;
    inputSet.insert(tgtGenome);
    inputSet.insert(mrca);
    getGenomesInSpanningTree(inputSet, pathSet);
    genomesOnPath = &pathSet;
  }

  hal_size_t numResults = DefaultMappedSegment::map(this, outSegments,
                                                    tgtGenome,
                                                    genomesOnPath, doDupes,
                                                    minLength,
                                                    coalescenceLimit,
                                                    mrca);
  return numResults;
}

void DefaultSegmentIterator::print(ostream& os) const
{
  hal_index_t ai = getArrayIndex();

  os << "gen=" << getGenome()->getName()
     << " seq=" << getSequence()->getName()
     << " idx=" << ai;

  bool offRight = 
     isTop() ? ai >= (hal_index_t)getGenome()->getNumTopSegments() :
     ai >= (hal_index_t)getGenome()->getNumBottomSegments();

  if (ai != NULL_INDEX && !offRight)
  {
    os 
       << " start=" << getStartPosition()
       << " end=" << getEndPosition()
       << " len=" << getLength()
       << " off=[" << getStartOffset() << "," << getEndOffset() << "]"
       << " rev=" << getReversed();
  }
}

//////////////////////////////////////////////////////////////////////////////
// SLICED SEGMENT INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultSegmentIterator::toReverse() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
}

void DefaultSegmentIterator::toReverseInPlace() const
{
  assert (inRange() == true);
  _reversed = !_reversed;
  swap(_startOffset, _endOffset);
}

hal_offset_t DefaultSegmentIterator::getStartOffset() const
{
  return _startOffset;
}

hal_offset_t DefaultSegmentIterator::getEndOffset() const
{
  return _endOffset;
}

void DefaultSegmentIterator::slice(hal_offset_t startOffset, 
                                   hal_offset_t endOffset) const
{
  assert(startOffset < getSegment()->getLength());
  assert(endOffset < getSegment()->getLength());
  _startOffset = startOffset;
  _endOffset = endOffset;
}


bool DefaultSegmentIterator::getReversed() const
{
  return _reversed;
}
   
//////////////////////////////////////////////////////////////////////////////
// SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void DefaultSegmentIterator::toLeft(hal_index_t leftCutoff) const
{
  if (_reversed == false)
  {
    if (_startOffset == 0)
    {
      getSegment()->setArrayIndex(getGenome(), 
                                  getSegment()->getArrayIndex() - 1);
      _endOffset = 0;
    }
    else
    {
      _endOffset = getSegment()->getLength() - _startOffset;
      _startOffset = 0;
    }
    if (getSegment()->getArrayIndex() >= 0 && 
        leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      assert(getSegment()->getStartPosition() <= leftCutoff);
      _startOffset = leftCutoff - getSegment()->getStartPosition();
    }
  }
  
  else
  {
    if (_startOffset == 0)
    {
      getSegment()->setArrayIndex(getGenome(), 
                                  getSegment()->getArrayIndex() + 1);
      _endOffset = 0;
    }
    else
    {
      _endOffset = getSegment()->getLength() - _startOffset;
      _startOffset = 0;
    }    
    if ((hal_size_t)getArrayIndex() < getNumSegmentsInGenome() &&
        leftCutoff != NULL_INDEX && overlaps(leftCutoff))
    {
      _startOffset = getSegment()->getStartPosition() + 
         getSegment()->getLength() - 1 - leftCutoff;
    }
  }
  assert((hal_size_t)getArrayIndex() >= 
         getNumSegmentsInGenome() ||
         getArrayIndex() < 0 || 
         _startOffset + _endOffset <= getSegment()->getLength());
}

void DefaultSegmentIterator::toRight(hal_index_t rightCutoff) const  
{
  if (_reversed == false)
  {
    if (_endOffset == 0)
    {
      getSegment()->setArrayIndex(getGenome(), 
                                  getSegment()->getArrayIndex() + 1);
      _startOffset = 0;
    }
    else
    {
      _startOffset =  getSegment()->getLength() - _endOffset;
      _endOffset = 0;
    }
    
    if ((hal_size_t)getArrayIndex() < getNumSegmentsInGenome() &&
        rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = getSegment()->getStartPosition() +
         getSegment()->getLength() - rightCutoff - 1;
    }
  }
  
  else
  {
    if (_endOffset == 0)
    {
      getSegment()->setArrayIndex(getGenome(), 
                                  getSegment()->getArrayIndex() - 1);
      _startOffset = 0;
    }
    else
    {
      _startOffset =  getSegment()->getLength() - _endOffset;
      _endOffset = 0;
    }

    if (rightCutoff != NULL_INDEX && overlaps(rightCutoff))
    {
      _endOffset = rightCutoff - getSegment()->getStartPosition(); 
    }
  }
  assert ((hal_size_t)getArrayIndex() >= 
          getNumSegmentsInGenome() ||
          getArrayIndex() < 0 || 
          _startOffset + _endOffset <= getSegment()->getLength());
}

void DefaultSegmentIterator::toSite(hal_index_t position, bool slice) const
{
  const Genome* genome = getGenome();
  hal_index_t len = (hal_index_t)genome->getSequenceLength();
  hal_index_t nseg = (hal_index_t)getNumSegmentsInGenome();
  
  assert(len != 0);
  double avgLen = (double)len / (double)nseg;
  hal_index_t hint = 
     (hal_index_t)min(nseg - 1., avgLen * ((double)position / (double)len));
  getSegment()->setArrayIndex(genome, hint);
  _startOffset = 0;
  _endOffset = 0;
  
  // out of range
  if (position < 0)
  {
    getSegment()->setArrayIndex(genome, NULL_INDEX);
    return;
  }
  else if (position >= len)
  {
    getSegment()->setArrayIndex(genome, len);
    return;
  }

  hal_index_t left = 0;
  hal_index_t leftStartPosition = 0;
  hal_index_t right = nseg - 1;
  hal_index_t rightStartPosition = len - 1;
  assert(getSegment()->getArrayIndex()  >= 0 &&  
         getSegment()->getArrayIndex() < nseg);

  while (overlaps(position) == false)
  {
    assert(left != right);
    if (rightOf(position) == true)
    {
      right = getSegment()->getArrayIndex();
      rightStartPosition = getSegment()->getStartPosition();
      avgLen = double(rightStartPosition - leftStartPosition) / (right - left);
      hal_index_t delta = (hal_index_t)
         max((rightStartPosition - position) / avgLen, 1.);
      delta = min(delta, getSegment()->getArrayIndex());
      getSegment()->setArrayIndex(genome, 
                                  getSegment()->getArrayIndex() - delta);
      assert(getSegment()->getArrayIndex()  >= 0 &&
             getSegment()->getArrayIndex() < nseg);
    }
    else
    {
      assert(leftOf(position) == true);
      left = getSegment()->getArrayIndex();
      leftStartPosition = getSegment()->getStartPosition();
      avgLen = double(rightStartPosition - leftStartPosition) / (right - left);
      hal_index_t delta = (hal_index_t)
         max((position - leftStartPosition) / avgLen, 1.);
      delta = min(delta, nseg - 1 - getSegment()->getArrayIndex());
      getSegment()->setArrayIndex(genome,
                                  getSegment()->getArrayIndex() + delta);
      assert(getSegment()->getArrayIndex() >= 0 &&
             getSegment()->getArrayIndex() < nseg);
    }
  }
  
  assert(overlaps(position) == true);
  
  if (slice == true)
  {
    _startOffset = position - getSegment()->getStartPosition();
    _endOffset = getSegment()->getStartPosition() + 
       getSegment()->getLength() - position - 1;
    if (_reversed)
    {
//       swap(_startOffset, _endOffset);
    }
  }  
}

