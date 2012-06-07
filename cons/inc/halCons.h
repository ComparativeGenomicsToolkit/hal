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

namespace hal {

struct ConsStats
{
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
};

class HalCons
{
public:

   HalCons();
   HalCons(hal::AlignmentConstPtr alignment); 
   virtual ~HalCons();

   void printCsv(std::ostream& outStream) const;
   void analyzeAlignment(hal::AlignmentConstPtr alignment);

protected:

   void analyzeGenomeRecursive(const std::string& genomeName);

   typedef std::pair<std::string, std::string> StrPair;
   typedef std::map<StrPair, ConsStats> BranchMap;

   BranchMap _branchMap;
   hal::AlignmentConstPtr _alignment;

};

}

inline std::ostream& operator<<(std::ostream& os, const hal::HalCons& halCons)
{
  halCons.printCsv(os);
  return os;
}

#endif
