/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halBlockLiftover.h"

using namespace std;
using namespace hal;

BlockLiftover::BlockLiftover() : Liftover()
{

}

BlockLiftover::~BlockLiftover()
{

}

void BlockLiftover::visitBegin()
{
  if (_srcGenome->getNumTopSegments() > 0)
  {
    _refSeg = _srcGenome->getTopSegmentIterator();
    _lastIndex = (hal_index_t)_srcGenome->getNumTopSegments();
  }
  else
  {
    _refSeg = _srcGenome->getBottomSegmentIterator();
    _lastIndex = (hal_index_t)_srcGenome->getNumBottomSegments();
  }

  set<const Genome*> inputSet;
  inputSet.insert(_srcGenome);
  inputSet.insert(_tgtGenome);
  getGenomesInSpanningTree(inputSet, _spanningTree);
}

void BlockLiftover::liftInterval()
{
  hal_index_t globalStart = _bedLine._start + _srcSequence->getStartPosition();
  hal_index_t globalEnd = _bedLine._end - 1 + _srcSequence->getStartPosition();
  
  _refSeg->toSite(globalStart, false);
  hal_offset_t startOffset = globalStart - _refSeg->getStartPosition();
  hal_offset_t endOffset = 0;
  if (globalEnd <= _refSeg->getEndPosition())
  {
    endOffset = _refSeg->getEndPosition() - globalEnd;
  }
  _refSeg->slice(startOffset, endOffset);
  
  assert(_refSeg->getStartPosition() ==  globalStart);
  assert(_refSeg->getEndPosition() <= globalEnd);

  while (_refSeg->getArrayIndex() < _lastIndex &&
         _refSeg->getStartPosition() <= globalEnd)  
  {
    if (_bedLine._strand == '-')
    {
      _refSeg->toReverseInPlace();
    }
    _refSeg->getMappedSegments(_mappedSegments, _tgtGenome, &_spanningTree,
                              _traverseDupes);
    if (_bedLine._strand == '-')
    {
      _refSeg->toReverseInPlace();
    }
    _refSeg->toRight(globalEnd);
  }

  // need to clean up set here

  for (std::set<MappedSegmentConstPtr>::iterator i = _mappedSegments.begin();
       i != _mappedSegments.end(); ++i)
  {
    const Sequence* seq = (*i)->getSequence();
    hal_size_t seqStart = seq->getStartPosition();
    _outBedLines.push_back(_bedLine);
    BedLine& outBedLine = _outBedLines.back();
    outBedLine._chrName = seq->getName();
    outBedLine._start = (*i)->getStartPosition();
    outBedLine._end = (*i)->getEndPosition() + 1 - seqStart;
    if ((*i)->getReversed() == true)
    {
      outBedLine._strand = '-';
      swap(outBedLine._start, outBedLine._end);
    }
    else
    {
      outBedLine._strand = '+';
    }
    assert(outBedLine._start < outBedLine._end);
  }
}

