/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halBlockLiftover.h"
#include "halBlockMapper.h"

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
    _refSeg->getMappedSegments(_mappedSegments, _tgtGenome, &_spanningTree,
                              _traverseDupes);
    _refSeg->toRight(globalEnd);
  }

  cleanTargetParalogies();
  vector<MappedSegmentConstPtr> fragments;
  BlockMapper::MSSet emptySet;
  set<hal_index_t> queryCutSet;
  set<hal_index_t> targetCutSet;
  
  for (std::set<MappedSegmentConstPtr>::iterator i = _mappedSegments.begin();
       i != _mappedSegments.end(); ++i)
  {
    BlockMapper::extractSegment(i, emptySet, fragments, &_mappedSegments, 
                                targetCutSet, queryCutSet);

    const Sequence* seq = (*i)->getSequence();
    hal_size_t seqStart = seq->getStartPosition();
    _outBedLines.push_back(_bedLine);
    BedLine& outBedLine = _outBedLines.back();
    outBedLine._chrName = seq->getName();
    outBedLine._start = min(min(fragments.front()->getStartPosition(), 
                                fragments.front()->getEndPosition()),
                            min(fragments.back()->getStartPosition(),
                                fragments.back()->getEndPosition()));
    outBedLine._start -= seqStart;
    outBedLine._end = max(max(fragments.front()->getStartPosition(), 
                                fragments.front()->getEndPosition()),
                            max(fragments.back()->getStartPosition(),
                                fragments.back()->getEndPosition()));
    outBedLine._end -= seqStart;
    outBedLine._strand = (*i)->getReversed() ? '-' : '+';
    assert(outBedLine._start < outBedLine._end);
  }
}

void BlockLiftover::cleanTargetParalogies()
{
  set<MappedSegmentConstPtr>::iterator i;
  set<MappedSegmentConstPtr>::iterator j;
  for (i = _mappedSegments.begin(); i != _mappedSegments.end(); i = j)
  {
    j = i;
    ++j;
    if (j != _mappedSegments.end())
    {
      if ((*i)->getStartPosition() == (*j)->getStartPosition() &&
          (*i)->getEndPosition() == (*j)->getEndPosition())
      {
        assert((*i)->getSequence() == (*j)->getSequence());
        _mappedSegments.erase(i);
      }
    }
  }
}

