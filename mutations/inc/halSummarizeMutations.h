/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMUTATIONS_H
#define _HALMUTATIONS_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"
#include "halAverage.h"
#include "halMutationsStats.h"

namespace hal {

class SummarizeMutations
{
public:

   SummarizeMutations();
   virtual ~SummarizeMutations();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(hal::AlignmentConstPtr alignment,
                         hal_size_t gapThreshold,
                         double nThreshold,
                         bool justSubs,
                         const std::set<std::string>* targetSet = NULL);

protected:

   void analyzeGenomeRecursive(const std::string& genomeName);   
   void substitutionAnalysis(const Genome* genome, MutationsStats& stats);
   void rearrangementAnalysis(const Genome* genome, MutationsStats& stats);
   void subsAndGapInserts(GappedTopSegmentIteratorConstPtr gappedTop, 
                          MutationsStats& stats);

   typedef std::pair<std::string, std::string> StrPair;
   typedef std::map<StrPair, MutationsStats> BranchMap;

   BranchMap _branchMap;
   hal::AlignmentConstPtr _alignment;
   hal_size_t _gapThreshold;
   double _nThreshold;
   bool _justSubs;
   const std::set<std::string>* _targetSet;
};

}

inline std::ostream& operator<<(std::ostream& os, 
                                const hal::SummarizeMutations& halCons)
{
  halCons.printCsv(os);
  return os;
}

#endif
