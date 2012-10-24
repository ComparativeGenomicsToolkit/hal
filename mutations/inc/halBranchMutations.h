/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBRANCHMUTATIONS_H
#define _HALBRANCHMUTATIONS_H

#include <iostream>
#include <string>
#include <map>
#include "hal.h"
#include "halMutationsStats.h"

namespace hal {

class BranchMutations
{
public:

   BranchMutations();
   virtual ~BranchMutations();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(AlignmentConstPtr alignment,
                         hal_size_t gapThreshold,
                         std::ostream* snpStream,
                         std::ostream* svStream,
                         const SegmentedSequence* reference,
                         hal_index_t startPosition,
                         hal_size_t length,
                         const std::set<const Genome*>* targets);

protected:

   typedef std::map<const Genome*, PositionCache*> CacheMap;
   typedef std::map<const Genome*, MutationsStats> StatsMap;

   AlignmentConstPtr _alignment;
   std::ostream* _snpStream;
   std::ostream* _svStream;
   hal_size_t _maxGap;
   ColumnIteratorConstPtr _colIt;
   CacheMap _cacheMap;
   StatsMap _statsMap;
};

}

#endif
