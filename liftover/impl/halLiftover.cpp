/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halLiftover.h"

using namespace std;
using namespace hal;

Liftover::Liftover() : _outBedStream(NULL),                       
                       _inBedVersion(-1), _outBedVersion(-1),
                       _srcGenome(NULL), _tgtGenome(NULL)
{

}

Liftover::~Liftover()
{

}

void Liftover::convert(AlignmentConstPtr alignment,
                       const Genome* srcGenome,
                       istream* inBedStream,
                       const Genome* tgtGenome,
                       ostream* outBedStream,
                       bool addExtraColumns,
                       bool traverseDupes)
{
  _srcGenome = srcGenome;
  _tgtGenome = tgtGenome;
  _outBedStream = outBedStream;
  _addExtraColumns = addExtraColumns;
  _traverseDupes = traverseDupes;
  _missedSet.clear();
  _tgtSet.clear();
  assert(_srcGenome && inBedStream && tgtGenome && outBedStream);

  _tgtSet.insert(tgtGenome);
  
  if (_inBedVersion == -1)
  {
    _inBedVersion = BedScanner::getBedVersion(inBedStream);
  }
  if (_outBedVersion == -1)
  {
    _outBedVersion = _inBedVersion;
  }

  scan(inBedStream, _inBedVersion);
}

void Liftover::visitLine()
{
  _outBedLines.clear();
  _srcSequence = _srcGenome->getSequence(_bedLine._chrName);
  if (_srcSequence == NULL)
  {
    pair<set<string>::iterator, bool> result = _missedSet.insert(
      _bedLine._chrName);
    if (result.second == true)
    {
      std::cerr << "Unable to find sequence " << _bedLine._chrName 
                << " in genome " << _srcGenome->getName() << endl;
    }
  }
      
  else if (_bedLine._end > (hal_index_t)_srcSequence->getSequenceLength())
  {
    std::cerr << "Skipping interval with endpoint " << _bedLine._end 
              << "because sequence " << _bedLine._chrName << " has length " 
              << _srcSequence->getSequenceLength() << endl;
  }
  
  if (_inBedVersion > 9 && !_bedLine._blocks.empty())
  {
    liftBlockIntervals();
  }
  else
  {
    liftInterval();
  }
  writeLineResults();
}

void Liftover::writeLineResults()
{
  if (_outBedVersion > 9)
  {
    collapseExtendedBedLines();
  }
  BedList::iterator i = _outBedLines.begin();
  for (; i != _outBedLines.end(); ++i)
  {
    if (_addExtraColumns == false)
    {
      i->_extra.clear();
    }
    i->write(*_outBedStream, _outBedVersion);
  }
}

void Liftover::collapseExtendedBedLines()
{
  _outBedLines.sort();
  
   // want to greedily merge up colinear intervals
  BedList::iterator i = _outBedLines.begin();
  BedList::iterator j;
  BedList::iterator k;
  
  for (; i != _outBedLines.end(); ++i)
  {
    j = i;
    ++j;
    for (; j != _outBedLines.end() && 
            i->_chrName == j->_chrName &&
            i->_strand == j->_strand; j = k)
    {
      k = j;
      k++;
      if (canMerge(*i, *j) == true)
      {
        mergeAsBlockInterval(*i, *j);
        _outBedLines.erase(j);
      }
    }
  }   
}

void Liftover::liftBlockIntervals()
{
  BedLine blockBed = _bedLine;
  std::sort(blockBed._blocks.begin(), blockBed._blocks.end());
  _bedLine._blocks.clear();
  hal_index_t prev = _bedLine._start;
  vector<BedBlock>::iterator blockIt = blockBed._blocks.begin();
  for (; blockIt != blockBed._blocks.end(); ++blockIt)
  {
    // region before this block
    _bedLine._start = prev;
    _bedLine._end = blockIt->_start;
    assert(_bedLine._end >= _bedLine._start);
    if (_bedLine._end > _bedLine._start)
    {
      liftInterval();
    }
    // the block
    _bedLine._start = blockIt->_start;
    _bedLine._end = blockIt->_start + blockIt->_length;
    if (_bedLine._end > _bedLine._start)
    {
      liftInterval();
    }
    prev = _bedLine._end;
  }
  
  // bit between last block and end
  _bedLine._start = prev;
  _bedLine._end = blockBed._end;
  if (_bedLine._start > _bedLine._end)
  {
    liftInterval();
  }
}

bool Liftover::canMerge(const BedLine& bedLine1, const BedLine& bedLine2)
{
  assert(bedLine2._blocks.empty());
  assert(bedLine1._chrName == bedLine2._chrName);
  assert(bedLine1._strand == bedLine2._strand);

  return (bedLine1._blocks.empty() ||
          bedLine2._start >= bedLine1._blocks.back()._start +
          bedLine1._blocks.back()._length);
}

void Liftover::mergeAsBlockInterval(BedLine& bedTarget, 
                                    const BedLine& bedSource)
{
  assert(canMerge(bedTarget, bedSource) == true);
  if (bedTarget._blocks.empty())
  {
    BedBlock block = {bedTarget._start, bedTarget._end - bedTarget._start};
    bedTarget._blocks.push_back(block);
  }

  BedBlock block = {bedSource._start, bedSource._end - bedSource._start};
  bedTarget._start = std::min(bedTarget._start, bedSource._start);
  bedTarget._end = std::max(bedTarget._end, bedSource._end);
  bedTarget._blocks.push_back(block);
}
