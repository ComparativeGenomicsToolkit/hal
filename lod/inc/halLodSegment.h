/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLODSEGMENT_H
#define _HALLODSEGMENT_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <list>
#include <map>
#include <cstdlib>
#include <cmath>
#include "hal.h"

namespace hal {

class LodSegment;
class LodBlock;

std::ostream& operator<<(std::ostream& os, const LodSegment& segment);

struct LodSegmentPLess
{
   bool operator()(const LodSegment* s1, const LodSegment* s2) const;
};

/* Segment strcture for the Level of Detail graph.  These Segments are what 
 * will eventually be turned into (top or bottom) segments in the output
 * HAL graph. Segments are oriented (have a head and a tail).  They have
 * a tail adjacency and a head adjacency to connect to other segments on
 * the same sequence.  Segments with head before tail are called flipped.
 */
class LodSegment
{
   friend std::ostream& operator<<(std::ostream& os,const LodSegment& segment);

public:

   LodSegment();
   LodSegment(LodBlock* block, const Sequence* sequence, 
              hal_index_t pos, bool flipped);
   ~LodSegment();

   // inline get methods
   const Sequence* getSequence() const;
   hal_index_t getTailPos() const;
   hal_index_t getHeadPos() const;
   hal_index_t getLeftPos() const;
   hal_index_t getRightPos() const;
   bool getFlipped() const;
   hal_size_t getLength() const;
   const LodSegment* getTailAdj() const;
   const LodSegment* getHeadAdj() const;
   hal_size_t getTailAdjLen() const;
   hal_size_t getHeadAdjLen() const;
   bool getTailToHead() const;
   bool getTailToTail() const;
   bool getHeadToTail() const;
   bool getHeadToHead() const;
   bool overlaps(const LodSegment& other) const;
   hal_index_t getArrayIndex() const;
   void setArrayIndex(hal_index_t index);
   LodBlock* getBlock() const;

   /** Add a new edge from right endpoint of this segment to left
    * endpoint of tgt segment.  (whether or not these endpoints are 
    * heads or tails depends on the orientation of the segments)
    */
   void addEdgeFromRightToLeft(LodSegment* tgt);

   /** Extend the segment tail-wards by extLen */
   void extendTail(hal_size_t extLen);

   /** Extend the segment head-wards by extLen */
   void extendHead(hal_size_t extLen);

   /** Insert a new segment as a head adjacency to this
    * segment. So the head of this will be connected to the 
    * tail of the new segment, and the new segment's head
    * will connect to whatever this's head connected to. The
    * new segment will have 0 distance from this segment, and
    * it's length is given by the parameter.  The new segment
    * is then returned */
   LodSegment* insertNewHeadAdj(LodBlock* block, hal_size_t newLen);
   LodSegment* insertNewTailAdj(LodBlock* block, hal_size_t newLen);

   /** Merge the head adjacency segment to this segment.  That segment
    * should then get taken out of consideration */
   void mergeHead();
   
protected:


   /* Find the minimum edge lengths in both directions */
   void getMinLengths(hal_size_t& rMin, hal_size_t& fMin) const;

   const Sequence* _sequence;
   hal_index_t _tailPos;
   hal_index_t _afterHeadPos;
   LodSegment* _tailAdj;
   LodSegment* _headAdj;
   hal_index_t _arrayIndex;
   LodBlock* _block;

private:
   LodSegment(const LodSegment&);
   const LodSegment& operator=(const LodSegment&) const;
};

inline bool LodSegmentPLess::operator()(const LodSegment* s1, 
                                        const LodSegment* s2) const
{
  bool left = s1->getRightPos() < s2->getLeftPos();
  return left;
}

inline const Sequence* LodSegment::getSequence() const
{
  return _sequence;
}

inline hal_index_t LodSegment::getTailPos() const
{
  return _tailPos;
}

inline hal_index_t LodSegment::getHeadPos() const
{
  return getFlipped() ? _afterHeadPos + 1 : _afterHeadPos - 1;
}

inline hal_index_t LodSegment::getLeftPos() const
{
  return getFlipped() ? getHeadPos() : getTailPos();
}

inline hal_index_t LodSegment::getRightPos() const
{
  return getFlipped() ? getTailPos() : getHeadPos();
}

inline bool LodSegment::getFlipped() const
{
  assert(_afterHeadPos != _tailPos);
  return _afterHeadPos < _tailPos;
}

inline hal_size_t LodSegment::getLength() const
{
  return (hal_size_t)std::abs(_afterHeadPos - _tailPos);
}

inline const LodSegment* LodSegment::getTailAdj() const
{
  return _tailAdj;
}

inline const LodSegment* LodSegment::getHeadAdj() const
{
  return _headAdj;
}

inline hal_size_t LodSegment::getTailAdjLen() const
{
  assert(_tailAdj != NULL);
  hal_index_t otherPos = 
     getTailToHead() ? _tailAdj->getHeadPos() : _tailAdj->getTailPos();
  assert(otherPos > getRightPos() || otherPos < getLeftPos());
  return (hal_size_t)std::abs(getTailPos() - otherPos) - 1;
}

inline hal_size_t LodSegment::getHeadAdjLen() const
{
  assert(_headAdj != NULL);
  hal_index_t otherPos = 
     getHeadToTail() ? _headAdj->getTailPos() : _headAdj->getHeadPos();
  assert(otherPos > getRightPos() || otherPos < getLeftPos());
  return (hal_size_t)std::abs(getHeadPos() - otherPos) - 1;
}

inline bool LodSegment::getTailToHead() const
{
  assert(_tailAdj != NULL);
  assert(this == _tailAdj->_headAdj || this == _tailAdj->_tailAdj);
  return this == _tailAdj->_headAdj;
}

inline bool LodSegment::getTailToTail() const
{
  return !getTailToHead();
}

inline bool LodSegment::getHeadToTail() const
{
  assert(_headAdj != NULL);
  assert(this == _headAdj->_headAdj || this == _headAdj->_tailAdj);
  return this == _headAdj->_tailAdj;
}

inline bool LodSegment::getHeadToHead() const
{
  return !getHeadToTail();
}

inline bool LodSegment::overlaps(const LodSegment& other) const
{
  return other.getRightPos() >= getLeftPos() && 
     other.getLeftPos() <= getRightPos();
}
  
inline hal_index_t LodSegment::getArrayIndex() const
{
  return _arrayIndex;
}

inline void LodSegment::setArrayIndex(hal_index_t index)
{
  _arrayIndex = index;
}

inline LodBlock* LodSegment::getBlock() const
{
  return _block;
}

}

#endif
