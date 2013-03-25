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
  hal_size_t minTail = getMinTailAdjLen();
  minTail = (hal_size_t)std::ceil(maxFrac * (double)minTail);
  
  for (LodBlock::SegmentIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    (*i)->extendTail(minTail);
    assert((*i)->getLength() == getLength());
  }

  hal_size_t minHead = getMinHeadAdjLen();
  minHead = (hal_size_t)std::ceil(maxFrac * (double)minHead);
  
  for (LodBlock::SegmentIterator i = _segments.begin();
       i != _segments.end(); ++i)
  {
    (*i)->extendHead(minHead);
    assert((*i)->getLength() == getLength());
  }
}


hal_size_t LodBlock::getMinTailAdjLen() const
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

hal_size_t LodBlock::getMinHeadAdjLen() const
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
