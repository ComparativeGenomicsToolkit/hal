/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMUTATIONSSTATS_H
#define _HALMUTATIONSSTATS_H

#include <iostream>
#include <string>
#include <vector>
#include "hal.h"
#include "halAverage.h"

namespace hal {

struct MutationsStats
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
   hal_size_t _matches;

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

   static void printHeader(std::ostream& os);
};

std::ostream& operator<<(std::ostream& os, const MutationsStats& stats);

MutationsStats& operator+=(MutationsStats& ms, const MutationsStats& other);
MutationsStats& operator/=(MutationsStats& ms, hal_size_t N);

}

#endif
