/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"
#include "halBlockMapper.h"

using namespace hal;
using namespace std;

hal_size_t BlockMapper::_maxAdjScan = 3;

BlockMapper::BlockMapper()
{

}

BlockMapper::~BlockMapper()
{
  erase();
}

void BlockMapper::erase()
{
  _segMap.clear();
  _flipSet.clear();
  _spanningTree.clear();
}

void BlockMapper::init(const Genome* refGenome, const Genome* queryGenome,
                       hal_index_t absRefFirst, hal_index_t absRefLast,
                       bool doDupes, hal_size_t minLength,
                       bool mapTargetAdjacencies)
{
  erase();
  _absRefFirst = absRefFirst;
  _absRefLast = absRefLast;
  _doDupes = doDupes;
  _minLength = minLength;
  _mapAdj = mapTargetAdjacencies;
  _refGenome = refGenome;
  _refSequence = refGenome->getSequenceBySite(_absRefFirst);
  assert(_refSequence == refGenome->getSequenceBySite(_absRefLast));
  _queryGenome = queryGenome;

  set<const Genome*> inputSet;
  inputSet.insert(_refGenome);
  inputSet.insert(_queryGenome);
  getGenomesInSpanningTree(inputSet, _spanningTree);
}

void BlockMapper::map()
{
  SegmentIteratorConstPtr refSeg;
  hal_index_t lastIndex;
  set<const Genome*>  inputSet;
  inputSet.insert(_refGenome);
  inputSet.insert(_queryGenome);
  if (getLowestCommonAncestor(inputSet) == _refGenome)
  {
    refSeg = _refGenome->getBottomSegmentIterator();
    lastIndex = _refGenome->getNumBottomSegments();
  }
  else
  {
    refSeg = _refGenome->getTopSegmentIterator();
    lastIndex = _refGenome->getNumTopSegments();
  }

  refSeg->toSite(_absRefFirst, false);
  hal_offset_t startOffset = _absRefFirst - refSeg->getStartPosition();
  hal_offset_t endOffset = 0;
  if (_absRefLast <= refSeg->getEndPosition())
  {
    endOffset = refSeg->getEndPosition() - _absRefLast;
  }
  refSeg->slice(startOffset, endOffset);
  
  assert(refSeg->getStartPosition() ==  _absRefFirst);
  assert(refSeg->getEndPosition() <= _absRefLast);

  while (refSeg->getArrayIndex() < lastIndex &&
         refSeg->getStartPosition() <= _absRefLast)  
  {
    refSeg->getMappedSegments(_segMap, _queryGenome, &_spanningTree,
                              _doDupes, _minLength);
    refSeg->toRight(_absRefLast);
  }

  if (_mapAdj)
  {
    // construct a set sorted on query instead of ref. 
    _flipSet.clear();
    _flipSet.insert(_segMap.begin(), _segMap.end());

    MSFlipSet::const_iterator i;
    for (i = _flipSet.begin(); i != _flipSet.end(); ++i)
    {
      mapAdjacencies(i);
    }
  }
}

void BlockMapper::extractReferenceParalogies(MSSet& outParalogies)
{
  MSFlipSet::iterator i = _flipSet.begin();
  MSFlipSet::iterator j = _flipSet.end();
  MSFlipSet::iterator k;
  hal_index_t iStart = NULL_INDEX;
  hal_index_t iEnd = NULL_INDEX;
  bool iIns = false;
  hal_index_t jStart;
  hal_index_t jEnd;
  while (i != _flipSet.end())
  {
    // we divide list sorted on query segments into equivalence
    // classes based on their query coordinates.  these classes
    // are orderdered by their reference coordinates.  we keep
    // the lowest entry in each class (by ref coordinates) and 
    // extract all the others into outParalgoies
    if (iStart == NULL_INDEX)
    {
      iIns = false;
      iStart = (*i)->getStartPosition();
      iEnd = (*i)->getEndPosition();
      if ((*i)->getReversed())
      {
        swap(iStart, iEnd);
      }
    }
    j = i;
    ++j;
    while (j != _flipSet.end())
    {
      jStart = (*j)->getStartPosition();
      jEnd = (*j)->getEndPosition();
      if ((*j)->getReversed())
      {
        swap(jStart, jEnd);
      }
      // note we have overlaps here (possibly) that will need to 
      // be clipped down the line
      if (iStart <= jStart && iEnd >= jStart)
      {
        if (!iIns)
        {
          iIns = true;
          outParalogies.insert(*i);
        }
        outParalogies.insert(*j);
        k = j;
        ++k;
        _segMap.erase(*j);
        _flipSet.erase(j);
        j = k;
      }
      else
      {
        iStart = NULL_INDEX;
        iIns = false;
        break;
      }
    }
    i = j;
  }
}

  
void BlockMapper::mapAdjacencies(MSFlipSet::const_iterator flipIt)
{
  assert(_flipSet.empty() == false && flipIt != _flipSet.end());
  assert(_segMap.find(*flipIt) != _segMap.end());
  MappedSegmentConstPtr mappedQuerySeg = *flipIt;
  hal_index_t maxIndex;
  hal_index_t minIndex;
  SegmentIteratorConstPtr queryIt = makeIterator(mappedQuerySeg, 
                                                 minIndex,
                                                 maxIndex);
  MSSet backResults;
  MSFlipSet::const_iterator flipNext = flipIt;
  if (queryIt->getReversed())
  {
    flipNext = flipNext == _flipSet.begin() ? _flipSet.end() : --flipNext;
  }
  else
  {
    ++flipNext;
  }

  hal_size_t iter = 0;
  queryIt->toRight();
  while (queryIt->getArrayIndex() >= minIndex &&
         queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan)
  {
    bool wasCut = false;
    if (flipNext != _flipSet.end())
    {
      // if cut returns nothing, then the region in question is covered
      // by flipNext (ie already mapped).
      wasCut = cutByNext(queryIt, *flipNext, !queryIt->getReversed());
    }
    if (wasCut == true)
    {
      break;
    }
    size_t backSize = backResults.size();
    assert(queryIt->getArrayIndex() >= 0);
    queryIt->getMappedSegments(backResults, _refGenome, &_spanningTree,
                               _doDupes, _minLength);
    // something was found, that's good enough.
    if (backResults.size() > backSize)
    {
      break;
    }
    queryIt->toRight();
    ++iter;
  }

  queryIt = makeIterator(mappedQuerySeg, 
                         minIndex,
                         maxIndex);

  MSFlipSet::const_iterator flipPrev = flipIt;
  if (queryIt->getReversed())
  {
    ++flipPrev;
  }
  else
  {
    flipPrev = flipPrev == _flipSet.begin() ? _flipSet.end() : --flipPrev;
  }
  iter = 0;
  queryIt->toLeft();
  while (queryIt->getArrayIndex() >= minIndex &&
         queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan)
  {
    bool wasCut = false;
    if (flipPrev != _flipSet.end())
    {
      // if cut returns nothing, then the region in question is covered
      // by flipPrev (ie already mapped).
      wasCut = cutByNext(queryIt, *flipPrev, queryIt->getReversed());
    }
    if (wasCut == true)
    {
      break;
    }
    size_t backSize = backResults.size();
    queryIt->getMappedSegments(backResults, _refGenome, &_spanningTree,
                               _doDupes, _minLength);
    // something was found, that's good enough.
    if (backResults.size() > backSize)
    {
      break;
    }
    queryIt->toLeft();
    ++iter;
  }

  // flip the results and copy back to our main set.
  for (MSSet::iterator i = backResults.begin(); i != backResults.end(); ++i)
  {
    MappedSegmentConstPtr mseg = *i;
    if (mseg->getSequence() == _refSequence)
    {
      mseg->flip();
      SlicedSegmentConstPtr refSeg = mseg->getSource();
      if (refSeg->getReversed())
      {
        mseg->fullReverse();
      }

      MSSet::const_iterator j = _flipSet.lower_bound(*i);
      bool overlaps = false;
      if (j != _flipSet.begin())
      {
        --j;
      }
      for (size_t count = 0; count < 3 && j != _flipSet.end() && !overlaps; 
           ++count, ++j)
      {
        overlaps = mseg->overlaps((*j)->getStartPosition()) ||
           mseg->overlaps((*j)->getEndPosition()) || 
           (*j)->overlaps(mseg->getStartPosition()) ||
           (*j)->overlaps(mseg->getEndPosition());
      }
      if (overlaps == false)
      {
        _segMap.insert(mseg);
      }
    }
  }
}

SegmentIteratorConstPtr BlockMapper::makeIterator(
  MappedSegmentConstPtr mappedSegment, hal_index_t& minIndex,
  hal_index_t& maxIndex)
{
  SegmentIteratorConstPtr segIt;
  if (mappedSegment->isTop())
  {
    segIt = mappedSegment->getGenome()->getTopSegmentIterator(
      mappedSegment->getArrayIndex());
    minIndex = segIt->getSequence()->getTopSegmentArrayIndex();
    maxIndex = minIndex + 
       (hal_index_t)segIt->getSequence()->getNumTopSegments();
  }
  else
  {
    segIt = mappedSegment->getGenome()->getBottomSegmentIterator(
      mappedSegment->getArrayIndex());
    minIndex = segIt->getSequence()->getBottomSegmentArrayIndex();
    maxIndex = minIndex + 
       (hal_index_t)segIt->getSequence()->getNumBottomSegments();
  }
  
  if (mappedSegment->getReversed())
  {
    segIt->toReverse();
  }
  segIt->slice(mappedSegment->getStartOffset(), 
               mappedSegment->getEndOffset());
  
  assert(segIt->getGenome() == mappedSegment->getGenome());
  assert(segIt->getArrayIndex() == mappedSegment->getArrayIndex());
  assert(segIt->getStartOffset() == mappedSegment->getStartOffset());
  assert(segIt->getEndOffset() == mappedSegment->getEndOffset());
  assert(segIt->getReversed() == mappedSegment->getReversed());

  return segIt;
}

bool BlockMapper::cutByNext(SlicedSegmentConstPtr queryIt, 
                            SlicedSegmentConstPtr nextSeg,
                            bool right)
{
  assert(queryIt->getGenome() == nextSeg->getGenome());
  assert(queryIt->isTop() == nextSeg->isTop());
  bool wasCut = false;

  if (queryIt->getArrayIndex() == nextSeg->getArrayIndex())
  {
    hal_offset_t so1 = queryIt->getStartOffset();
    hal_offset_t eo1 = queryIt->getEndOffset();
    if (queryIt->getReversed())
    {
      swap(so1, eo1);
    }
    hal_offset_t so2 = nextSeg->getReversed() ? nextSeg->getEndOffset() :
       nextSeg->getStartOffset();

    if (right)
    {
      // overlap on start position.  we zap
      if (so1 >= so2)
      {
        wasCut = true;
      }
      else
      {
        hal_index_t e1 = max(queryIt->getEndPosition(), 
                             queryIt->getStartPosition());
        hal_index_t s2 = min(nextSeg->getEndPosition(), 
                             nextSeg->getStartPosition());
        // end position of queryIt overlaps next seg.  so we cut it.
        if (e1 >= s2)
        {
          hal_index_t delta = 1 + e1 - s2;
          assert(delta < (hal_index_t)queryIt->getLength());
          hal_offset_t newEndOffset = eo1 + delta;
          hal_offset_t newStartOffset = so1;
          if (queryIt->getReversed() == true)
          {
            swap(newEndOffset, newStartOffset);
          }
          queryIt->slice(newStartOffset, newEndOffset);
        }
      }
    }
    else
    {
      hal_index_t s1 = min(queryIt->getEndPosition(), 
                           queryIt->getStartPosition());
      hal_index_t e1 = max(queryIt->getEndPosition(), 
                           queryIt->getStartPosition());
      hal_index_t e2 = max(nextSeg->getEndPosition(), 
                           nextSeg->getStartPosition());

      // overlap on end position.  we zap
      if (e1 <= e2)
      {
        wasCut = true;
      }
      else
      {
        // end position of queryIt overlaps next seg.  so we cut it.
        if (s1 <= e2)
        {
          hal_index_t delta = 1 + e2 - s1;
          assert(delta < (hal_index_t)queryIt->getLength());
          hal_offset_t newStartOffset = so1 + delta;
          hal_offset_t newEndOffset = eo1;
          if (queryIt->getReversed() == true)
          {
            swap(newEndOffset, newStartOffset);
          }
          queryIt->slice(newStartOffset, newEndOffset);
        }
      }
    }
  }
  return wasCut;
}

void BlockMapper::extractSegment(MSSet::iterator start, 
                                 const MSSet& paraSet,
                                 vector<MappedSegmentConstPtr>& fragments,
                                 MSSet* startSet)
{
  fragments.clear();
  fragments.push_back(*start);
  MSSet::iterator next = start;
  const Sequence* startSeq = (*start)->getSequence();
  MSSet::iterator temp;
  MappedSegmentConstPtr prev = *start;
  ++next;
  MSSet::key_compare comp = startSet->key_comp();
  bool keepOn = true;

  while (keepOn)
  {
    keepOn = false;
    SlicedSegmentConstPtr pseg = prev->getSource();
    assert(pseg->getReversed() == false);
    assert(next == startSet->end() || 
           (*next)->getSource()->getReversed() == false);

    while (next != startSet->end() && 
           (*next)->getSource()->getStartPosition() <=
           pseg->getEndPosition() + 1)
    {
      assert((*next)->getSource()->getReversed() == false);
      bool para1 = paraSet.find(prev) != paraSet.end();
      bool para2 = paraSet.find(*next) != paraSet.end();
      if (para1 == para2 &&
          (*next)->getSource()->getStartPosition() == 
          pseg->getEndPosition() + 1 && 
          (*next)->getSequence() == startSeq &&
          canMergeBlock(prev, *next) == true)
      {
        fragments.push_back(*next);
        prev = *next;
        temp = next;
        ++next;
        startSet->erase(temp);
        keepOn = true;
        break;
      }
      prev = *next;
      ++next;
    }
  }
  assert(fragments.front()->getSequence() == fragments.back()->getSequence());
}

bool BlockMapper::canMergeBlock(MappedSegmentConstPtr firstQuerySeg, 
                                MappedSegmentConstPtr lastQuerySeg)
{
  //return false;
  SlicedSegmentConstPtr firstRefSeg = firstQuerySeg->getSource();
  SlicedSegmentConstPtr lastRefSeg = lastQuerySeg->getSource();
  assert(firstRefSeg->getReversed() == false);
  assert(lastRefSeg->getReversed() == false);
  assert(firstRefSeg->getSequence() == lastRefSeg->getSequence());
  assert(firstQuerySeg->getGenome() == lastQuerySeg->getGenome());

  if (lastQuerySeg->getStartPosition() == _absRefFirst ||
      firstQuerySeg->getEndPosition() == _absRefLast)
  {
    return false;
  }
  if (firstQuerySeg->getReversed() == lastQuerySeg->getReversed() &&
      lastRefSeg->getStartPosition() == firstRefSeg->getEndPosition() + 1 &&
      firstQuerySeg->getSequence() == lastQuerySeg->getSequence())
  {
    if (firstQuerySeg->getReversed() == false &&
        lastQuerySeg->getStartPosition() == 
        firstQuerySeg->getEndPosition() + 1)
    {
      return true;
    }
    else if (firstQuerySeg->getReversed() == true &&
             lastQuerySeg->getStartPosition() == 
             firstQuerySeg->getEndPosition() - 1)
    {
      return true;
    }
  }
  return false;
}
