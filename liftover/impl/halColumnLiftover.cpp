/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halColumnLiftover.h"

using namespace std;
using namespace hal;

ColumnLiftover::ColumnLiftover() : Liftover()
{

}

ColumnLiftover::~ColumnLiftover()
{

}


void ColumnLiftover::liftInterval()
{  
  PositionMap posCacheMap;
  _colIt = _srcSequence->getColumnIterator(&_tgtSet, 0, _bedLine._start, 
                                           _bedLine._end - 1,
                                           !_traverseDupes);
  while (true) 
  {
    const ColumnMap* cMap = _colIt->getColumnMap();
    for (ColumnMap::const_iterator i = cMap->begin(); i != cMap->end(); ++i)
    {
      if (i->first->getGenome() == _tgtGenome)
      {
        const DNASet* dSet = i->second;
        const Sequence* seq = i->first;
        // if we're not adding the column, don't bother keeping track
        SeqIndex seqIdx(seq, 0);
        for (DNASet::const_iterator j = dSet->begin(); j != dSet->end(); ++j)
        {
          pair<PositionMap::iterator, bool> res =
             posCacheMap.insert(pair<SeqIndex, PositionCache*>(seqIdx, NULL));
          if (res.second == true)
          {
            res.first->second = new PositionCache();
          }
          res.first->second->insert((*j)->getArrayIndex());
        }
      }
    }
    if (_colIt->lastColumn() == true)
    {
      break;
    }
    _colIt->toRight();
  } 

  PositionMap::iterator pcmIt;
  for (pcmIt = posCacheMap.begin(); pcmIt != posCacheMap.end(); ++pcmIt)
  {
    const Sequence* seq = pcmIt->first.first;
    _outParalogy = pcmIt->first.second;
    hal_size_t seqStart = seq->getStartPosition();
    PositionCache* posCache = pcmIt->second;
    const IntervalSet* iSet = posCache->getIntervalSet();
    for (IntervalSet::const_iterator k = iSet->begin(); k != iSet->end(); ++k)
    {
      _outBedLines.push_back(_bedLine);
      BedLine& outBedLine = _outBedLines.back();
      outBedLine._chrName = seq->getName();
      outBedLine._start = k->second - seqStart;
      outBedLine._end = k->first + 1 - seqStart;
    }
    delete posCache;
  }
}

