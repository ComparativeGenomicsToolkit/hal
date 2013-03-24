/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cmath>
#include <cassert>
#include "halLodSegment.h"
#include "halLodEdge.h"

using namespace std;
using namespace hal;

LodSegment::LodSegment() : _sequence(NULL), _tailPos(NULL_INDEX), 
                           _afterHeadPos(NULL_INDEX), _tailAdj(NULL),
                           _headAdj(NULL)
{

}

LodSegment::LodSegment(const Sequence* sequence, hal_index_t pos, 
                       bool flipped) :
  _sequence(sequence), 
  _tailPos(pos),
  _afterHeadPos(pos),
  _tailAdj(NULL),
  _headAdj(NULL)
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

void LodSegment::addEdgeFromRightToLeft(const Sequence* sequence, 
                                        LodSegment* tgt)
{
  assert(tgt != NULL);
  LodSegment*& myCap = getFlipped() ? _tailAdj : _headAdj;
  LodSegment*& tgtCap = tgt->getFlipped() ? tgt->_headAdj : tgt->_tailAdj;
  assert(myCap == NULL);
  assert(tgtCap == NULL);
  myCap = tgt;
  tgtCap = this;
}

ostream& hal::operator<<(ostream& os, const LodSegment& segment)
{
  os << "seg: " << &segment << " "
     << segment.getSequence()->getName()
     << " [" << segment.getLeftPos() 
     << " (" << (segment.getFlipped() ? "H" : "T") << ")"
     << ", " << segment.getRightPos()
     << " (" << (segment.getFlipped() ? "T" : "H") << ")"
     << "]"
     << " tAdj: " << segment._tailAdj << " len=" << segment.getTailAdjLen()
     << "(" << (segment.getTailToHead() ? "H" : "T")
     << " hAdj: " << segment._headAdj << " len=" << segment.getHeadAjdLen()
     << "(" << (segment.getHeadToTail() ? "T" : "H");
  return os;
}
