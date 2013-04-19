/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cmath>
#include <cassert>
#include <limits>
#include "halLodBlock.h"

using namespace std;
using namespace hal;

LodBlock::LodBlock()
{

}

LodBlock::~LodBlock()
{
  clear();
}

void LodBlock::addSegment(LodSegment* segment)
{
  assert(segment != NULL);
  assert(_segments.empty() || getLength() == segment->getLength());
  _segments.push_back(segment);
}

void LodBlock::clear()
{
  for (SegmentIterator i = _segments.begin(); i != _segments.end(); ++i)
  {
    delete *i;
  }
  _segments.clear();
}

hal_size_t LodBlock::getTotalAdjLength() const
{
  hal_size_t total = 0;
  for (LodBlock::SegmentConstIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    if ((*i)->getTailAdj() != NULL)
    {
      total += (*i)->getTailAdjLen();
    }
    if ((*i)->getHeadAdj() != NULL)
    {
      total += (*i)->getHeadAdjLen();
    }
  }
  return total;
}

void LodBlock::extend(double maxFrac)
{
  hal_size_t tailExtLen = getMaxTailExtensionLen();
  tailExtLen = (hal_size_t)std::ceil(maxFrac * (double)tailExtLen);
  
  for (LodBlock::SegmentIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    (*i)->extendTail(tailExtLen);
    assert((*i)->getLength() == getLength());
  }

  hal_size_t headExtLen = getMaxHeadExtensionLen();
  headExtLen = (hal_size_t)std::ceil(maxFrac * (double)headExtLen);
  
  for (LodBlock::SegmentIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    (*i)->extendHead(headExtLen);
    assert((*i)->getLength() == getLength());
  }
}

LodBlock* LodBlock::getHeadMergePartner()
{
  LodBlock* adjBlock = NULL;
  bool uniqueValidAdj = true;
  for (LodBlock::SegmentIterator i = _segments.begin();
       i != _segments.end() && uniqueValidAdj == true; ++i)
  {
    const LodSegment* headAdj = (*i)->getHeadAdj();

    // Dont' want to merge from or a telomere
    // or to a telomere
    // or if there is a non-zero adjacency 
    // or if there is a head to head adjacency
    if (headAdj == NULL || 
        headAdj->getHeadAdj() == NULL ||
        (*i)->getHeadAdjLen() != 0 ||
        (*i)->getHeadToTail() == false)
    {
      uniqueValidAdj = false;
    }
    else if (adjBlock == NULL && 
             headAdj->getBlock() != this &&
             headAdj->getBlock()->getNumSegments() == getNumSegments())
    {
      adjBlock = headAdj->getBlock();
    }
    else
    {
      uniqueValidAdj = adjBlock == headAdj->getBlock();
    }
  }
  return uniqueValidAdj ? adjBlock : NULL;
}

void LodBlock::mergeHead(LodBlock* adjBlock)
{
  assert(adjBlock == getHeadMergePartner());
  if (adjBlock != NULL)
  {
    for (LodBlock::SegmentIterator i = _segments.begin();
         i != _segments.end(); ++i)
    {
      (*i)->mergeHead();
    }
    adjBlock->clear();
  }
}

void LodBlock::insertNeighbours(vector<LodBlock*>& outList)
{
  LodBlock* lodBlock = NULL;
  while (true) 
  {
    lodBlock = insertNewTailNeighbour();
    if (lodBlock != NULL)
    {
      outList.push_back(lodBlock);
    }
    else
    {
      break;
    }
  }
  while (true) 
  {
    lodBlock = insertNewHeadNeighbour();
    if (lodBlock != NULL)
    {
      outList.push_back(lodBlock);
    }
    else
    {
      break;
    }
  }
  assert(getMaxHeadInsertionLen() == 0);
  assert(getMaxTailInsertionLen() == 0);
}

LodBlock* LodBlock::insertNewTailNeighbour()
{
  LodBlock* newBlock = NULL;
  hal_size_t maxTailInsLen = getMaxTailInsertionLen();
  if (maxTailInsLen > 0)
  {
    newBlock = new LodBlock();
    for (LodBlock::SegmentIterator i = _segments.begin();
         i != _segments.end(); ++i)
    {
      if ((*i)->getTailAdjLen()  >= maxTailInsLen)
      {
        LodSegment* newSeg = (*i)->insertNewTailAdj(newBlock, maxTailInsLen);
        newBlock->_segments.push_back(newSeg);
      }
    }
  }
  return newBlock;
}

LodBlock* LodBlock::insertNewHeadNeighbour()
{
  LodBlock* newBlock = NULL;
  hal_size_t maxHeadInsLen = getMaxHeadInsertionLen();
  if (maxHeadInsLen > 0)
  {
    newBlock = new LodBlock();
    for (LodBlock::SegmentIterator i = _segments.begin();
         i != _segments.end(); ++i)
    {
      if ((*i)->getHeadAdjLen() >= maxHeadInsLen)
      {
        LodSegment* newSeg = (*i)->insertNewHeadAdj(newBlock, maxHeadInsLen);
        newBlock->_segments.push_back(newSeg);
      }
    }
  }
  return newBlock;
}

hal_size_t LodBlock::getMaxTailExtensionLen() const
{
  hal_size_t minTail = numeric_limits<hal_size_t>::max();

  set<const LodSegment*> lookup;
  hal_size_t adjLen;
  for (LodBlock::SegmentConstIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    adjLen = (*i)->getTailAdjLen();
    if ((*i)->getTailToTail() && 
        lookup.find((*i)->getTailAdj()) != lookup.end())
    {
      // we have a tail to tail edge in the block.  can only extend half
      adjLen /= 2;
    }
    minTail = min(adjLen, minTail);
    lookup.insert(*i);
  }
  return minTail;
}

hal_size_t LodBlock::getMaxHeadExtensionLen() const
{
  hal_size_t minHead = numeric_limits<hal_size_t>::max();

  set<const LodSegment*> lookup;
  hal_size_t adjLen;
  for (LodBlock::SegmentConstIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    adjLen = (*i)->getHeadAdjLen();
    if ((*i)->getHeadToHead() && 
        lookup.find((*i)->getHeadAdj()) != lookup.end())
    {
      // we have a head to head edge in the block.  can only extend half
      adjLen /= 2;
    }
    minHead = min(adjLen, minHead);
    lookup.insert(*i);
  }
  return minHead;
}

hal_size_t LodBlock::getMaxTailInsertionLen() const
{
  hal_size_t minNZTail = numeric_limits<hal_size_t>::max();
  hal_size_t adjLen;
  for (LodBlock::SegmentConstIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    adjLen = (*i)->getTailAdjLen();
    if (adjLen > 0)
    {
      minNZTail = min(adjLen, minNZTail);
    }
  }
  return minNZTail == numeric_limits<hal_size_t>::max() ? 0 : minNZTail;
}

hal_size_t LodBlock::getMaxHeadInsertionLen() const
{
  hal_size_t minNZHead = numeric_limits<hal_size_t>::max();
  hal_size_t adjLen;
  for (LodBlock::SegmentConstIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    adjLen = (*i)->getHeadAdjLen();
    if (adjLen > 0)
    {
      minNZHead = min(adjLen, minNZHead);
    }
  }
  return minNZHead == numeric_limits<hal_size_t>::max() ? 0 : minNZHead;
}

ostream& hal::operator<<(ostream& os, const LodBlock& block)
{
  os << "block " << &block << " n=" << block.getNumSegments() 
     << " l=" << block.getLength() << ":";
  for (LodBlock::SegmentConstIterator i = block._segments.begin();
       i != block._segments.end(); ++i)
  {
    os << "\n  " << **i;
  }
  return os;
}
