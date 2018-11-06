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
#include "halMappedSegmentContainers.h"

namespace hal{

// helper class to make blocks for snake display
// keeps a map of segments in reference genome to sets of segments
// in query genome that map to it.  
// can optionally toggle dupes
class BlockMapper
{
public:

   BlockMapper();
   virtual ~BlockMapper();

   void init(const Genome* refGenome, const Genome* queryGenome, 
             hal_index_t absRefFirst, hal_index_t absRefLast,
             bool targetReversed,
             bool doDupes, hal_size_t minLength,
             bool mapTargetAdjacencies,
             const Genome *coalescenceLimit = NULL);
   void map();
   void extractReferenceParalogies(MappedSegmentSet& outParalogies);

   const MappedSegmentSet& getMap() const;
   MappedSegmentSet& getMap();

   static void extractSegment(MappedSegmentSet::iterator start, 
                              const MappedSegmentSet& paraSet,                       
                              std::vector<MappedSegmentPtr>& fragments,
                              MappedSegmentSet* startSet,
                              const std::set<hal_index_t>& targetCutPoints,
                              std::set<hal_index_t>& queryCutPoints);

   hal_index_t getAbsRefFirst() const;
   hal_index_t getAbsRefLast() const;

protected:
   
   void erase();
   void mapAdjacencies(MappedSegmentSet::const_iterator setIt);

   static SegmentIteratorPtr makeIterator(
     MappedSegmentPtr mappedSegment, 
     hal_index_t& minIndex,
     hal_index_t& maxIndex);

   // note queryIt passed by reference.  the pointer can be modified. 
   // ick.
   static bool cutByNext(SlicedSegmentConstPtr queryIt, 
                         SlicedSegmentConstPtr nextSeg,
                         bool right);
      
   static bool equalTargetStart(const MappedSegmentPtr& s1,
                                const MappedSegmentPtr& s2);
protected:

   MappedSegmentSet _segSet;
   MappedSegmentSet _adjSet;
   std::set<const Genome*> _downwardPath;
   std::set<const Genome*> _upwardPath;
   const Genome* _refGenome;
   const Sequence* _refSequence;
   const Genome* _queryGenome;
   hal_index_t _refChildIndex;
   hal_index_t _queryChildIndex;
   bool _doDupes;
   hal_size_t _minLength;
   hal_index_t _absRefFirst;
   hal_index_t _absRefLast;
   bool _mapAdj;
   bool _targetReversed;
   const Genome *_mrca;
   const Genome *_coalescenceLimit;

   static hal_size_t _maxAdjScan;
};

inline const MappedSegmentSet& BlockMapper::getMap() const
{
  return _segSet;
}

inline MappedSegmentSet& BlockMapper::getMap()
{
  return _segSet;
}

inline bool BlockMapper::equalTargetStart(const MappedSegmentPtr& s1,
                                          const MappedSegmentPtr& s2)
{
  hal_index_t p1 = std::min(s1->getStartPosition(), s1->getEndPosition());
  hal_index_t p2 = std::min(s2->getStartPosition(), s2->getEndPosition());
  assert(p1 != p2 ||
         std::max(s1->getStartPosition(), s1->getEndPosition()) ==
         std::max(s2->getStartPosition(), s2->getEndPosition()));
  return p1 == p2;
}

inline hal_index_t BlockMapper::getAbsRefFirst() const
{
  return _absRefFirst;
}

inline hal_index_t BlockMapper::getAbsRefLast() const
{
  return _absRefLast;
}

}

#endif
