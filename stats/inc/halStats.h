/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSTATS_H
#define _HALSTATS_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"

struct GenomeStats : public hal::Sequence::Info 
{
   size_t _numChildren;
   size_t _numSequences;
};

class HalStats
{
public:

   HalStats();
   HalStats(hal::AlignmentConstPtr alignment); 
   virtual ~HalStats();

   void printCsv(std::ostream& outStream) const;
   void readAlignment(hal::AlignmentConstPtr alignment);

protected:

   void readGenomeRecursive(const hal::Genome* genome);

   std::string _tree;
   std::vector<GenomeStats> _genomeStatsVec;
};

inline std::ostream& operator<<(std::ostream& os, const HalStats& halStats)
{
  halStats.printCsv(os);
  return os;
}

#endif
