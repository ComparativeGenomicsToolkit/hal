/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cmath>
#include <cassert>
#include "halLodSegment.h"

using namespace std;
using namespace hal;

LodSegment::LodSegment() : _sequence(NULL), _tailPos(NULL_INDEX), 
                           _afterHeadPos(NULL_INDEX), _tailAdj(NULL),
                           _headAdj(NULL), _arrayIndex(NULL_INDEX),
                           _block(NULL)
{

}

LodSegment::LodSegment(LodBlock* block, const Sequence* sequence, 
                       hal_index_t pos, bool flipped) :
  _sequence(sequence), 
  _tailPos(pos),
  _afterHeadPos(pos),
  _tailAdj(NULL),
  _headAdj(NULL),
  _arrayIndex(NULL_INDEX),
  _block(block)
{
  _afterHeadPos += flipped ? -1 : 1;
  assert(_sequence != NULL);
  assert(getLeftPos() >= getRightPos());
  assert(getLeftPos() >= _sequence->getStartPosition() - 1);
  assert(getRightPos() <= _sequence->getEndPosition() + 1);
}

LodSegment::~LodSegment()
{
}

void LodSegment::addEdgeFromRightToLeft(LodSegment* tgt)
{
  assert(tgt != NULL);
  LodSegment*& myCap = getFlipped() ? _tailAdj : _headAdj;
  LodSegment*& tgtCap = tgt->getFlipped() ? tgt->_headAdj : tgt->_tailAdj;
  assert(myCap == NULL);
  assert(tgtCap == NULL);
  myCap = tgt;
  tgtCap = this;
  assert(_tailAdj != _headAdj);
}

void LodSegment::extendTail(hal_size_t extLen)
{
  assert(getTailAdjLen() >= extLen);
  assert(getLeftPos() >= _sequence->getStartPosition() &&
         getRightPos() <= _sequence->getEndPosition());
  if (getFlipped() == true)
  {
    _tailPos += extLen;
  }
  else
  {
    _tailPos -= extLen;
  }
  assert(getLeftPos() >= _sequence->getStartPosition() &&
         getRightPos() <= _sequence->getEndPosition());
  assert(overlaps(*_tailAdj) == false);
}

void LodSegment::extendHead(hal_size_t extLen)
{
  assert(getHeadAdjLen() >= extLen);
  assert(getLeftPos() >= _sequence->getStartPosition() &&
         getRightPos() <= _sequence->getEndPosition());
  if (getFlipped() == true)
  {
    _afterHeadPos -= extLen;
  }
  else
  {
    _afterHeadPos += extLen;
  }
  assert(getLeftPos() >= _sequence->getStartPosition() &&
         getRightPos() <= _sequence->getEndPosition());
  assert(overlaps(*_headAdj) == false);
}

LodSegment* LodSegment::insertNewHeadAdj(LodBlock* block, hal_size_t newLen)
{
  assert(newLen > 0);
  hal_index_t newTailPos = getHeadPos();
  newTailPos += getFlipped() ? -1 : 1;
  LodSegment* newSeg = new LodSegment(block, getSequence(), newTailPos,
                                      getFlipped());
  bool headToHead = getHeadToHead();
  newSeg->_headAdj = _headAdj;
  if (headToHead)
  {
    assert(_headAdj->_headAdj == this);
    _headAdj->_headAdj = newSeg;
  }
  else
  {
    assert(_headAdj->_tailAdj == this);
    _headAdj->_tailAdj = newSeg;
  }
  newSeg->_tailAdj = this;
  _headAdj = newSeg;
  assert(getHeadAdjLen() == 0);
  assert(newSeg->getTailAdjLen() == 0);
  // we do the minus 1 below because the segment already
  // counts for 1
  newSeg->extendHead(newLen - 1);
  return newSeg;
}

LodSegment* LodSegment::insertNewTailAdj(LodBlock* block, hal_size_t newLen)
{
  assert(newLen > 0);
  hal_index_t newHeadPos = getTailPos();
  newHeadPos += getFlipped() ? 1 : -1;
  LodSegment* newSeg = new LodSegment(block, getSequence(), newHeadPos, 
                                      getFlipped());
  bool tailToTail = getTailToTail();
  newSeg->_tailAdj = _tailAdj;
  if (tailToTail)
  {
    assert(_tailAdj->_tailAdj == this);
    _tailAdj->_tailAdj = newSeg;
  }
  else
  {
    assert(_tailAdj->_headAdj == this);
    _tailAdj->_headAdj = newSeg;
  }
  newSeg->_headAdj = this;
  _tailAdj = newSeg;
  assert(getTailAdjLen() == 0);
  assert(newSeg->getHeadAdjLen() == 0);
  // we do the minus 1 below because the segment already
  // counts for 1
  newSeg->extendTail(newLen - 1);
  return newSeg;
}

void LodSegment::mergeHead()
{
  assert(getHeadAdjLen() == 0 && getHeadToTail() == true);
  assert(getArrayIndex() == NULL_INDEX);
  LodSegment* segToMerge = _headAdj;
  LodSegment* newHeadAdj = segToMerge->_headAdj;
  assert(segToMerge->getTailAdj() == this);
  bool newHeadToTail = segToMerge->getHeadToTail();
  hal_size_t segToMergeHeadAdjLen = segToMerge->getHeadAdjLen();
  (void)segToMergeHeadAdjLen;
  hal_size_t segToMergeLen = segToMerge->getLength();
  _headAdj = newHeadAdj;
  if (newHeadToTail)
  {
    newHeadAdj->_tailAdj = this;
  }
  else
  {
    newHeadAdj->_headAdj = this;
  }
  assert(getHeadAdjLen() == segToMergeLen + segToMergeHeadAdjLen);
  extendHead(segToMergeLen);
  assert(getHeadAdjLen() == segToMergeHeadAdjLen);
  segToMerge->_headAdj = NULL;
  segToMerge->_tailAdj = NULL;
}

ostream& hal::operator<<(ostream& os, const LodSegment& segment)
{
  os << "seg: " << &segment << " "
     << segment.getSequence()->getName()
     << " [" << segment.getLeftPos() 
     << "(" << (segment.getFlipped() ? "H" : "T") << ")"
     << ", " << segment.getRightPos()
     << "(" << (segment.getFlipped() ? "T" : "H") << ")"
     << "]"
     << "[ai=" << segment.getArrayIndex() << "]"
     << " tAdj: " << segment._tailAdj;
  if (segment._tailAdj != NULL)
  {
    os  << " len=" << segment.getTailAdjLen() 
        << "(" << (segment.getTailToHead() ? "H)" : "T)");
  }
  os << " hAdj: " << segment._headAdj;
  if (segment._headAdj != NULL)
  {
    os  << " len=" << segment.getHeadAdjLen()
        << "(" << (segment.getHeadToTail() ? "T)" : "H)");
  }
  return os;
}
