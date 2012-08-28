/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCONS_H
#define _HALCONS_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"
#include "halAverage.h"

namespace hal {

struct ConsStats
{
   // old stats
   hal_size_t _genomeLength;
   hal_size_t _parentLength;
   double _branchLength;
   hal_size_t _subs;
   hal_size_t _numInserts;
   hal_size_t _numInsertBases;
   hal_size_t _numDeletes;
   hal_size_t _numDeleteBases;
   hal_size_t _numInverts;
   hal_size_t _numInvertBases;
   hal_size_t _numDups;
   hal_size_t _numDupBases;

   // stats from new rearragnemtn 
   typedef Average<hal_size_t> Avg;
   Avg _nothingLength;
   Avg _inversionLength;
   Avg _insertionLength;
   Avg _deletionLength;
   Avg _transpositionLength;
   Avg _duplicationLength;
   Avg _otherLength;
   Avg _gapInsertionLength;
   Avg _gapDeletionLength;

};

class HalCons
{
public:

   HalCons(hal::AlignmentConstPtr alignment, hal_size_t gapThreshold); 
   virtual ~HalCons();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(hal::AlignmentConstPtr alignment);

protected:

   void analyzeGenomeRecursive(const std::string& genomeName);   
   void rearrangementAnalysis(const Genome* genome, ConsStats& stats);

   typedef std::pair<std::string, std::string> StrPair;
   typedef std::map<StrPair, ConsStats> BranchMap;

   BranchMap _branchMap;
   hal::AlignmentConstPtr _alignment;
   hal_size_t _gapThreshold;

};

}

inline std::ostream& operator<<(std::ostream& os, const hal::HalCons& halCons)
{
  halCons.printCsv(os);
  return os;
}

#endif
