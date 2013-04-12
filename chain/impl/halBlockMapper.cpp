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

hal_size_t BlockMapper::_maxAdjScan = 5;

BlockMapper::BlockMapper()
{

}

BlockMapper::~BlockMapper()
{
  erase();
}

void BlockMapper::erase()
{
  for (SegMap::iterator i = _segMap.begin(); i != _segMap.end(); ++i)
  {
    delete i->second;
  }
  _segMap.clear();
}

void BlockMapper::init(const Genome* refGenome, const Genome* queryGenome,
                       hal_index_t absRefFirst, hal_index_t absRefLast,
                       bool doDupes)
{
  erase();
  _absRefFirst = absRefFirst;
  _absRefLast = absRefLast;
  _doDupes = doDupes;
  _refGenome = refGenome;
  _refSequence = refGenome->getSequenceBySite(_absRefFirst);
  assert(_refSequence == refGenome->getSequenceBySite(_absRefLast));
  _queryGenome = queryGenome;

  if (queryGenome->getParent() == refGenome)
  {
    _rel = RefParent;
    _queryChildIndex = refGenome->getChildIndex(queryGenome);
    _refChildIndex = NULL_INDEX;
  }
  else if (refGenome->getParent() == queryGenome)
  {
    _rel = RefChild;
    _queryChildIndex = NULL_INDEX;
    _refChildIndex = queryGenome->getChildIndex(refGenome);
  }
  else if (queryGenome->getParent() == refGenome->getParent())
  {
    _rel = RefSister;
    const Genome* parentGenome = queryGenome->getParent();
    _queryChildIndex = parentGenome->getChildIndex(queryGenome);
    _refChildIndex = parentGenome->getChildIndex(refGenome);
  }
  else
  {
    stringstream ss;
    ss << "Query species " << queryGenome->getName() << " is neither the "
       << "parent, child, nor sibling of the Reference species " 
       << refGenome->getName();
    throw hal_exception(ss.str());
  }
}

void BlockMapper::map()
{
  SegmentIteratorConstPtr refSeg;
  hal_index_t lastIndex;
  if (_rel == RefParent)
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
    if (_segMap.find(refSeg) == _segMap.end())
    {
      mapRef(refSeg);
    }
    refSeg->toRight(_absRefLast);
  }

  // Copy all the reference segments we found because 
  // adding adjacencies will modify the structure
  vector<SegMap::iterator> firstRoundResults(_segMap.size());
  hal_size_t i = 0;
  for (SegMap::iterator mapIt = _segMap.begin(); mapIt != _segMap.end(); 
       ++mapIt)
  {
    firstRoundResults[i++] = mapIt;
  }
  // For every query for every referencec, update the map with adjacent
  // blocks when found
  for (i = 0; i < firstRoundResults.size(); ++i)
  {
    SegMap::iterator mapIt = firstRoundResults[i];
    for (SegSet::iterator setIt = mapIt->second->begin(); 
         setIt != mapIt->second->end(); ++setIt)
    {
      mapAdjacencies(*setIt);
    }
  }
}

void BlockMapper::addParalogies(TopSegmentIteratorConstPtr top, 
                                SegSet* segSet)
{
  if (_doDupes)
  {
    while (top->hasNextParalogy() == true &&
           top->getNextParalogyIndex() != (*segSet->begin())->getArrayIndex())
    {
      top->toNextParalogy();
      segSet->insert(top->copy());
    }
  }
}

bool BlockMapper::isCanonical(TopSegmentIteratorConstPtr top)
{
  BottomSegmentIteratorConstPtr bottom = 
     top->getGenome()->getParent()->getBottomSegmentIterator(
       top->getParentIndex());
  hal_index_t childIndex = top->getGenome() == _refGenome ? _refChildIndex :
     _queryChildIndex;
  return bottom->getChildIndex(childIndex) == top->getArrayIndex();
}

void BlockMapper::mapRef(SegmentIteratorConstPtr refSeg)
{
  assert(!refSeg->getReversed());
  if (_rel == RefParent)
  { 
    mapRefParent(refSeg);
  }
  else if (_rel == RefChild)
  {
      mapRefChild(refSeg);
  }
  else // _rel == RefSister
  {
    mapRefSister(refSeg);
  }
}

void BlockMapper::mapRefParent(SegmentIteratorConstPtr refSeg)
{
  BottomSegmentIteratorConstPtr refBottom = 
     refSeg.downCast<BottomSegmentIteratorConstPtr>();

  if (refBottom->hasChild(_queryChildIndex))
  {
    TopSegmentIteratorConstPtr queryTop = _queryGenome->getTopSegmentIterator();
    queryTop->toChild(refBottom, _queryChildIndex);
    SegSet* segSet = new SegSet();
    segSet->insert(queryTop->copy());
    addParalogies(queryTop, segSet);
    assert(_segMap.find(refBottom) == _segMap.end());
    _segMap.insert(pair<SegmentIteratorConstPtr, SegSet*>(
                     refBottom->copy(), segSet));
  }
}

void BlockMapper::mapRefChild(SegmentIteratorConstPtr refSeg)
{
  TopSegmentIteratorConstPtr refTop = 
     refSeg.downCast<TopSegmentIteratorConstPtr>();

  if (refTop->hasParent() && (_doDupes || isCanonical(refTop)))
  {
    BottomSegmentIteratorConstPtr queryBottom = 
       _queryGenome->getBottomSegmentIterator();
    queryBottom->toParent(refTop);
    SegSet* segSet = new SegSet();
    segSet->insert(queryBottom->copy());
    assert(_segMap.find(refTop) == _segMap.end());
    _segMap.insert(pair<SegmentIteratorConstPtr, SegSet*>(
                     refTop->copy(), segSet));
  }
}

void BlockMapper::mapRefSister(SegmentIteratorConstPtr refSeg)
{
  TopSegmentIteratorConstPtr refTop = 
     refSeg.downCast<TopSegmentIteratorConstPtr>();

  if (refTop->hasParent() && (_doDupes || isCanonical(refTop)))
  {
    BottomSegmentIteratorConstPtr parBottom = 
       _refGenome->getParent()->getBottomSegmentIterator();
    parBottom->toParent(refTop);
    if (parBottom->hasChild(_queryChildIndex))
    {
      TopSegmentIteratorConstPtr queryTop = 
         _queryGenome->getTopSegmentIterator();
      queryTop->toChild(parBottom, _queryChildIndex);
      SegSet* segSet = new SegSet();
      segSet->insert(queryTop->copy());
      addParalogies(queryTop, segSet);
      assert(_segMap.find(queryTop) == _segMap.end());
      _segMap.insert(pair<SegmentIteratorConstPtr, SegSet*>(
                       refTop->copy(), segSet));
    }
  }
}
  
void BlockMapper::mapAdjacencies(SegmentIteratorConstPtr querySeg)
{
  SegmentIteratorConstPtr querySegCpy;
  SegmentIteratorConstPtr querySegCpy2;

  hal_index_t maxIndex;
  hal_index_t minIndex;
  if (_rel == RefParent || _rel == RefSister)
  {
    TopSegmentIteratorConstPtr queryTop =
       querySeg.downCast<TopSegmentIteratorConstPtr>();
    querySegCpy = queryTop->copy();
    querySegCpy2 = queryTop->copy();
    minIndex = queryTop->getSequence()->getTopSegmentArrayIndex();
    maxIndex = minIndex + 
       (hal_index_t)queryTop->getSequence()->getNumTopSegments();
  }
  else
  {
    BottomSegmentIteratorConstPtr queryBottom =
       querySeg.downCast<BottomSegmentIteratorConstPtr>();
    querySegCpy = queryBottom->copy();
    querySegCpy2 = queryBottom->copy();
    minIndex = queryBottom->getSequence()->getBottomSegmentArrayIndex();
    maxIndex = minIndex + 
       (hal_index_t)queryBottom->getSequence()->getNumBottomSegments();
  }

  hal_size_t iter = 0;
  querySegCpy->toRight();
  while ((querySegCpy->getReversed() || 
          querySegCpy->getArrayIndex() < maxIndex) &&
         (!querySegCpy->getReversed() || 
          querySegCpy->getArrayIndex() >= minIndex) && 
         iter < _maxAdjScan)
  {
    SegmentIteratorConstPtr refSeg = getAdjacencyInRef(querySegCpy);
    // never want an inverted reference.  
    if (refSeg->getReversed())
    {
      refSeg->toReverse();
      //but we don't want to keep in forward coordinates
      refSeg->slice(refSeg->getEndOffset(), refSeg->getStartOffset());
    }

    if (refSeg.get() != NULL)
    {
      if (_segMap.find(refSeg) == _segMap.end())
      {
        mapRef(refSeg);
      }
      break;
    }
    querySegCpy->toRight();
    ++iter;
  }
  // typo prevention
  querySegCpy = SegmentIteratorConstPtr();

  iter = 0;
  querySegCpy2->toLeft();
  while ((querySegCpy2->getReversed() || 
          querySegCpy2->getArrayIndex() >= minIndex) &&
         (!querySegCpy2->getReversed() || 
          querySegCpy2->getArrayIndex() < maxIndex) && 
         iter < _maxAdjScan)
  {
    SegmentIteratorConstPtr refSeg = getAdjacencyInRef(querySegCpy2);
    // never want an inverted reference.  
    if (refSeg->getReversed())
    {
      refSeg->toReverse();
      //but we don't want to keep in forward coordinates
      refSeg->slice(refSeg->getEndOffset(), refSeg->getStartOffset());
    }

    if (refSeg.get() != NULL)
    {
      if (_segMap.find(refSeg) == _segMap.end())
      {
        mapRef(refSeg);
      }
      break;
    }
    querySegCpy2->toLeft();
    ++iter;
  }
}

SegmentIteratorConstPtr BlockMapper::getAdjacencyInRef(
  SegmentIteratorConstPtr querySeg)
{
  SegmentIteratorConstPtr refSeg;
  if (_rel == RefParent)
  { 
    TopSegmentIteratorConstPtr queryTop = 
       querySeg.downCast<TopSegmentIteratorConstPtr>();
    if (queryTop->hasParent())
    {
      BottomSegmentIteratorConstPtr refBottom = 
         _refGenome->getBottomSegmentIterator();
      refBottom->toParent(queryTop);
      refSeg = refBottom;
    }
  }
  else if (_rel == RefChild)
  {
    BottomSegmentIteratorConstPtr queryBottom = 
       querySeg.downCast<BottomSegmentIteratorConstPtr>();
    if (queryBottom->hasChild(_refChildIndex))
    {
      TopSegmentIteratorConstPtr refTop = 
         _refGenome->getTopSegmentIterator();
      refTop->toChild(queryBottom, _refChildIndex);
      refSeg = refTop;
    }
  }
  else // _rel == QuerySister
  {
    TopSegmentIteratorConstPtr queryTop = 
       querySeg.downCast<TopSegmentIteratorConstPtr>();
    if (queryTop->hasParent())
    {
      BottomSegmentIteratorConstPtr parBot = 
         queryTop->getGenome()->getParent()->getBottomSegmentIterator();
      parBot->toParent(queryTop);
      if (parBot->hasChild(_refChildIndex))
      {
        TopSegmentIteratorConstPtr refTop = 
         _refGenome->getTopSegmentIterator();
        refTop->toChild(parBot, _refChildIndex);
        refSeg = refTop;
      }
    }
  }
  if (refSeg.get() != NULL)
  {
    assert(refSeg->getGenome() == _refGenome);
    if (refSeg->getSequence() != _refSequence)
    {
      refSeg = SegmentIteratorConstPtr();
    }
  }
  return refSeg;
}
