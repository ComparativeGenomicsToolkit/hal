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

namespace hal {

struct ConsStats
{
   typedef Average<hal_size_t> Avg;
   // Tree Information
   hal_size_t _genomeLength;
   hal_size_t _parentLength;
   double _branchLength;

   // Subsitution Information
   hal_size_t _subs;
   hal_size_t _transitions;
   hal_size_t _transversions;

   // Rearrangement Information
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

class Mutations
{
public:

   Mutations(hal::AlignmentConstPtr alignment, hal_size_t gapThreshold); 
   virtual ~Mutations();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(hal::AlignmentConstPtr alignment);

protected:

   void analyzeGenomeRecursive(const std::string& genomeName);   
   void rearrangementAnalysis(const Genome* genome, ConsStats& stats);
   void subsAndGapInserts(GappedTopSegmentIteratorConstPtr gappedTop, 
                          ConsStats& stats);

   typedef std::pair<std::string, std::string> StrPair;
   typedef std::map<StrPair, ConsStats> BranchMap;

   BranchMap _branchMap;
   hal::AlignmentConstPtr _alignment;
   hal_size_t _gapThreshold;

};

}

inline std::ostream& operator<<(std::ostream& os, const hal::Mutations& halCons)
{
  halCons.printCsv(os);
  return os;
}

#endif
