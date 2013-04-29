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

using namespace std;
using namespace hal;

DefaultSegmentIterator::DefaultSegmentIterator(hal_offset_t startOffset, 
                                               hal_offset_t endOffset,
                                               bool reversed) :
  DefaultSlicedSegment(startOffset,
                       endOffset,
                       reversed)
{

}

DefaultSegmentIterator::~DefaultSegmentIterator()
{

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

