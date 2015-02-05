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
  _mrca = getLowestCommonAncestor(inputSet);
  if (_coalescenceLimit == NULL) {
      _coalescenceLimit = _mrca;
  }

  inputSet.clear();
  inputSet.insert(_coalescenceLimit);
  inputSet.insert(_tgtGenome);
  getGenomesInSpanningTree(inputSet, _downwardPath);
}

void BlockLiftover::liftInterval(BedList& mappedBedLines)
{
  _mappedSegments.clear();
  hal_index_t globalStart = _bedLine._start + _srcSequence->getStartPosition();
  hal_index_t globalEnd = _bedLine._end - 1 + _srcSequence->getStartPosition();
  bool flip = _bedLine._strand == '-';

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
    if (flip == true)
    {
      _refSeg->toReverseInPlace();
    }
    _refSeg->getMappedSegments(_mappedSegments, _tgtGenome, &_downwardPath,
                               _traverseDupes, 0, _coalescenceLimit, _mrca);
    if (flip == true)
    {
      _refSeg->toReverseInPlace();
    }
    _refSeg->toRight(globalEnd);
  }

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
    mappedBedLines.push_back(_bedLine);
    BedLine& outBedLine = mappedBedLines.back();
    outBedLine._blocks.clear();
    outBedLine._chrName = seq->getName();
    outBedLine._start = min(min(fragments.front()->getStartPosition(), 
                                fragments.front()->getEndPosition()),
                            min(fragments.back()->getStartPosition(),
                                fragments.back()->getEndPosition()));
    outBedLine._start -= seqStart;
    outBedLine._end = 1 + max(max(fragments.front()->getStartPosition(), 
                                  fragments.front()->getEndPosition()),
                              max(fragments.back()->getStartPosition(),
                                  fragments.back()->getEndPosition()));
    outBedLine._end -= seqStart;
    outBedLine._strand = (*i)->getReversed() ? '-' : '+';

    SlicedSegmentConstPtr srcFront = fragments.front()->getSource();
    SlicedSegmentConstPtr srcBack = fragments.back()->getSource();
    outBedLine._srcStart = min(min(srcFront->getStartPosition(), 
                                   srcFront->getEndPosition()),
                               min(srcBack->getStartPosition(),
                                   srcBack->getEndPosition()));
    outBedLine._srcStrand = srcFront->getReversed() ? '-' : '+';

    if (_bedLine._strand == '.')
    {
      outBedLine._strand = '.';
      outBedLine._srcStrand = '.';
    }

    assert(outBedLine._start < outBedLine._end);

    if (_outPSL == true && !fragments.empty())
    {
      readPSLInfo(fragments, outBedLine);
    }
  }
}

void BlockLiftover::readPSLInfo(vector<MappedSegmentConstPtr>& fragments, 
                                BedLine& outBedLine)
{
  const Sequence* srcSequence = fragments[0]->getSource()->getSequence();
  const Sequence* tSequence = fragments[0]->getSequence();

  outBedLine._psl.resize(1);
  PSLInfo& psl = outBedLine._psl[0];
  psl._matches = 0;
  psl._misMatches = 0;
  psl._repMatches = 0;
  psl._nCount = 0;
  psl._qNumInsert = 0;
  psl._qBaseInsert = 0;
  psl._tNumInsert = 0;
  psl._tBaseInsert = 0;
  psl._qSeqName = srcSequence->getName();
  psl._qSeqSize = srcSequence->getSequenceLength();
  psl._qStrand = fragments[0]->getSource()->getReversed() ? '-' : '+';
  assert(outBedLine._srcStart >= srcSequence->getStartPosition());
  psl._qChromOffset = srcSequence->getStartPosition();
  psl._qEnd = outBedLine._srcStart + 
     (outBedLine._end - outBedLine._start);
  psl._tSeqSize = tSequence->getSequenceLength();
  psl._qBlockStarts.clear();

  string sBuf;
  string tBuf;
  for (size_t i = 0; i < fragments.size(); ++i)
  {
    assert(fragments[i]->getSource()->getReversed() == 
           fragments[0]->getSource()->getReversed());
    assert(fragments[i]->getSource()->getSequence() == srcSequence);
    assert(fragments[i]->getSequence() == tSequence);

    fragments[i]->getSource()->getString(sBuf);
    fragments[i]->getString(tBuf);
    
    for (size_t j = 0; j < sBuf.length(); ++j)
    {
      if (sBuf[j] == tBuf[j])
      {
        if (!isMasked(sBuf[j]) && !isMasked(tBuf[j]))
        {
          ++psl._matches;
        }
        else
        {
          ++psl._repMatches;
        }
      }
      else
      {
        ++psl._misMatches;
      }
      if (isMissingData(tBuf[j]))
      {
        ++psl._nCount;
      }
    }
  }
}
