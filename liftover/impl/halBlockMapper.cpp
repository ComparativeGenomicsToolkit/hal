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

hal_size_t BlockMapper::_maxAdjScan = 1;

BlockMapper::BlockMapper()
{

}

BlockMapper::~BlockMapper()
{
  erase();
}

void BlockMapper::erase()
{
  _segSet.clear();
  _adjSet.clear();
  _downwardPath.clear();
  _upwardPath.clear();
}

void BlockMapper::init(const Genome* refGenome, const Genome* queryGenome,
                       hal_index_t absRefFirst, hal_index_t absRefLast,
                       bool targetReversed,
                       bool doDupes, hal_size_t minLength,
                       bool mapTargetAdjacencies,
                       const Genome* coalescenceLimit)
{
  erase();
  _absRefFirst = absRefFirst;
  _absRefLast = absRefLast;
  _targetReversed = targetReversed;
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
  _mrca = getLowestCommonAncestor(inputSet);

  if (coalescenceLimit == NULL)
  {
    _coalescenceLimit = _mrca;
  }
  else
  {
    _coalescenceLimit = coalescenceLimit;
  }

  // The path between the coalescence limit (the highest point in the
  // tree) and the query genome is needed to traverse down into the
  // correct children.
  inputSet.clear();
  inputSet.insert(_queryGenome);
  inputSet.insert(_coalescenceLimit);
  getGenomesInSpanningTree(inputSet, _downwardPath);

  // similarly, the upward path is needed to get the adjacencies properly.
  inputSet.clear();
  inputSet.insert(_refGenome);
  inputSet.insert(_coalescenceLimit);
  getGenomesInSpanningTree(inputSet, _upwardPath);
}

void BlockMapper::map()
{
  SegmentIteratorConstPtr refSeg;
  hal_index_t lastIndex;
  if ((_mrca == _refGenome) && (_refGenome != _queryGenome))
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
    if (_targetReversed == true)
    {
      refSeg->toReverseInPlace();
    }
    refSeg->getMappedSegments(_segSet, _queryGenome, &_downwardPath,
                              _doDupes, _minLength, _coalescenceLimit, _mrca);
    if (_targetReversed == true)
    {
      refSeg->toReverseInPlace();            
    }
    refSeg->toRight(_absRefLast);
  }

  if (_mapAdj)
  {
    assert(_targetReversed == false);
    MSSet::const_iterator i;
    for (i = _segSet.begin(); i != _segSet.end(); ++i)
    {
      if (_adjSet.find(*i) == _adjSet.end())
      {
        mapAdjacencies(i);
      }
    }
  }
}

void BlockMapper::extractReferenceParalogies(MSSet& outParalogies)
{
  MSSet::iterator i = _segSet.begin();
  MSSet::iterator j = _segSet.end();
  MSSet::iterator k;
  hal_index_t iStart = NULL_INDEX;
  hal_index_t iEnd = NULL_INDEX;
  bool iIns = false;
  hal_index_t jStart;
  hal_index_t jEnd;
  while (i != _segSet.end())
  {
    // we divide list sorted on query segments into equivalence
    // classes based on their query coordinates.  these classes
    // are ordered by their reference coordinates.  we keep
    // the lowest entry in each class (by ref coordinates) and 
    // extract all the others into outParalogies
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
    // if (_adjSet.find(*i) == _adjSet.end())
    {
      while (j != _segSet.end())
      {
        jStart = (*j)->getStartPosition();
        jEnd = (*j)->getEndPosition();
        if ((*j)->getReversed())
        {
          swap(jStart, jEnd);
        }
        // note we should not have overlaps here because the mappedSegment
        // results cuts everything to be same or disjoint
        if (iStart == jStart)
        {
          assert(iEnd == jEnd);
          if (!iIns)
          {
            iIns = true;
            outParalogies.insert(*i);
          }
          outParalogies.insert(*j);
          k = j;
          ++k;
          _segSet.erase(*j);
          j = k;
        }
        else
        {
          assert(iStart > jEnd || iEnd < jStart);
          iStart = NULL_INDEX;
          iIns = false;
          break;
        }
      }
    }
    i = j;
  }

#ifndef _NDEBUG
  for (i = _segSet.begin(); i != _segSet.end(); ++i)
  {
    j = i;
    ++j;
    if (j!= _segSet.end())
    {
      iEnd = max((*i)->getEndPosition(), (*i)->getStartPosition());
      jStart = min((*j)->getStartPosition(), (*j)->getEndPosition());
      assert(jStart > iEnd);
    }
  }
#endif
}

  
void BlockMapper::mapAdjacencies(MSSet::const_iterator segIt)
{
  assert(_segSet.empty() == false && segIt != _segSet.end());
  MappedSegmentConstPtr mappedQuerySeg = *segIt;
  hal_index_t maxIndex;
  hal_index_t minIndex;
  SegmentIteratorConstPtr queryIt = makeIterator(mappedQuerySeg, 
                                                 minIndex,
                                                 maxIndex);
  MSSet backResults;
  MSSet::const_iterator segNext = segIt;
  if (queryIt->getReversed())
  {
    segNext = segNext == _segSet.begin() ? _segSet.end() : --segNext;
  }
  else
  {
    ++segNext;
  }

  hal_size_t iter = 0;
  queryIt->toRight();
  while (queryIt->getArrayIndex() >= minIndex &&
         queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan)
  {
    bool wasCut = false;
    if (segNext != _segSet.end())
    {
      // if cut returns nothing, then the region in question is covered
      // by segNext (ie already mapped).
      wasCut = cutByNext(queryIt, *segNext, !queryIt->getReversed());
    }
    if (wasCut == true)
    {
      break;
    }
    size_t backSize = backResults.size();
    assert(queryIt->getArrayIndex() >= 0);
    queryIt->getMappedSegments(backResults, _refGenome, &_upwardPath,
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

  MSSet::const_iterator segPrev = segIt;
  if (queryIt->getReversed())
  {
    ++segPrev;
  }
  else
  {
    segPrev = segPrev == _segSet.begin() ? _segSet.end() : --segPrev;
  }
  iter = 0;
  queryIt->toLeft();
  while (queryIt->getArrayIndex() >= minIndex &&
         queryIt->getArrayIndex() < maxIndex && iter < _maxAdjScan)
  {
    bool wasCut = false;
    if (segPrev != _segSet.end())
    {
      // if cut returns nothing, then the region in question is covered
      // by segPrev (ie already mapped).
      wasCut = cutByNext(queryIt, *segPrev, queryIt->getReversed());
    }
    if (wasCut == true)
    {
      break;
    }
    size_t backSize = backResults.size();
    queryIt->getMappedSegments(backResults, _refGenome, &_upwardPath,
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

      MSSet::const_iterator j = _segSet.lower_bound(*i);
      bool overlaps = false;
      if (j != _segSet.begin())
      {
        --j;
      }
      for (size_t count = 0; count < 3 && j != _segSet.end() && !overlaps; 
           ++count, ++j)
      {
        overlaps = mseg->overlaps((*j)->getStartPosition()) ||
           mseg->overlaps((*j)->getEndPosition()) || 
           (*j)->overlaps(mseg->getStartPosition()) ||
           (*j)->overlaps(mseg->getEndPosition());
      }
      if (overlaps == false)
      {
        _segSet.insert(mseg);
        _adjSet.insert(mseg);
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
                                 MSSet* startSet,
                                 const set<hal_index_t>& targetCutPoints,
                                 set<hal_index_t>& queryCutPoints)
{
  fragments.clear();
  fragments.push_back(*start);
  const Sequence* startSeq = (*start)->getSequence();

  vector<MSSet::iterator> vector1;
  vector<MSSet::iterator>* v1 = &vector1;
  vector<MSSet::iterator> vector2;
  vector<MSSet::iterator>* v2 = &vector2;
  vector<MSSet::iterator> toErase;

  v1->push_back(start);
  MSSet::iterator next = start;
  ++next;

  // equivalence class based on start set
  while (next != startSet->end() && equalTargetStart(*v1->back(), *next))
  {
    assert((*next)->getLength() == (*v1->back())->getLength());
    v1->push_back(next);
    ++next;
  }
  
  while (next != startSet->end())
  {
    // equivalence class based on next element
    while (next != startSet->end() && 
           (v2->empty() || equalTargetStart(*v2->back(), *next)) &&
           v2->size() < v1->size())
    {
      assert(v2->empty() || (*next)->getLength() == (*v2->back())->getLength());
      v2->push_back(next);
      ++next;
    }
    
    // check if all elements of each class are compatible for a merge
    assert(v1->size() > 0 && v2->size() > 0);
    bool canMerge = v1->size() == v2->size();
    for (size_t i = 0; i < v1->size() && canMerge; ++i)
    {    
      canMerge = 
         (v1->size() == v2->size() && (*v2->at(i))->getSequence() == startSeq &&
          (*v1->at(i))->canMergeRightWith(*v2->at(i), &queryCutPoints,
                                          &targetCutPoints) && 
          (paraSet.find(*v1->at(i)) == paraSet.end()) == 
          (paraSet.find(*v2->at(i)) == paraSet.end()));
    }
    if (canMerge == true)
    {
      // add the next element and flag for deletion from start set
      fragments.push_back(*v2->at(0));
      toErase.push_back(v2->at(0));
    }
    else
    {
      break;
    }
    v1->clear();
    swap(v1, v2);
  }
  
  // we extracted a paralgous region along query coordinates.  we 
  // add its right query coordinate to the cutset to make sure none
  // of the other paralogs ever get merged beyond it.
  if (v1->size() > 1)
  {
    queryCutPoints.insert(max(fragments.back()->getStartPosition(),
                              fragments.back()->getEndPosition()));
  }
   
  for (size_t i = 0; i < toErase.size(); ++i)
  {
    startSet->erase(toErase[i]);
  }

  assert(fragments.front()->getSequence() == fragments.back()->getSequence());
}

