/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBLOCKMAPPER_H
#define _HALBLOCKMAPPER_H

#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <map>
#include "hal.h"

namespace hal{

struct SegmentPtrLess
{
   bool operator()(const SegmentIteratorConstPtr& s1, 
                   const SegmentIteratorConstPtr& s2) const;
};

// helper class to make blocks for snake display
// keeps a map of segments in reference genome to sets of segments
// in query genome that map to it.  
// can optionally toggle dupes
class BlockMapper
{
public:

   typedef std::set<SegmentIteratorConstPtr, SegmentPtrLess> SegSet;
   typedef std::map<SegmentIteratorConstPtr, SegSet*, SegmentPtrLess> SegMap;

   BlockMapper();
   virtual ~BlockMapper();

   void init(const Genome* refGenome, const Genome* queryGenome, 
             hal_index_t absRefFirst, hal_index_t absRefLast,
             bool doDupes);
   void map();

   const SegMap& getMap() const;

protected:
   
   enum Relation {RefParent, RefChild, RefSister};

   void erase();
   void addParalogies(TopSegmentIteratorConstPtr top, SegSet* segSet);
   bool isCanonical(TopSegmentIteratorConstPtr top);

   void mapRef(SegmentIteratorConstPtr refSeg);
   void mapRefParent(SegmentIteratorConstPtr refSeg);
   void mapRefChild(SegmentIteratorConstPtr refSeg);
   void mapRefSister(SegmentIteratorConstPtr refSeg);

   void mapAdjacencies(SegmentIteratorConstPtr querySeg);
   SegmentIteratorConstPtr getAdjacencyInRef(
     SegmentIteratorConstPtr querySeg);
   
protected:

   SegMap _segMap;
   Relation _rel;
   const Genome* _refGenome;
   const Sequence* _refSequence;
   const Genome* _queryGenome;
   hal_index_t _refChildIndex;
   hal_index_t _queryChildIndex;
   bool _doDupes;
   hal_index_t _absRefFirst;
   hal_index_t _absRefLast;

   static hal_size_t _maxAdjScan;
};

inline bool
SegmentPtrLess::operator()(const SegmentIteratorConstPtr& s1, 
                           const SegmentIteratorConstPtr& s2) const
{
  return s1->getStartPosition() < s2->getStartPosition();
}

inline const BlockMapper::SegMap& BlockMapper::getMap() const
{
  return _segMap;
}

}

#endif
