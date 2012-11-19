/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halMutationsStats.h"

using namespace std;
using namespace hal;

void MutationsStats::printHeader(std::ostream& os)
{
  os << "BranchLength, GenomeLength," 
     " ParentLength, Subtitutions, Transitions, Transversions, Matches"
     " GapInsertions, GapInsertedBases, GapDeletions, GapDeletedBases," 
     " Insertions, InsertionBases, Deletions, DeletionBases, Inversions,"
     " InvertedBases, Duplications, DuplicatedBases, Transpositions,"
     " TranspositionBases, Other";
}

ostream& hal::operator<<(ostream& os, const MutationsStats& stats)
{
  // right now the break pairs of all detected events get marked as
  // other, so we subtract them from the total as they were already
  // counted.  (this may result in a slight underestimation of other
  // in some cases which will need to be looked into later
  hal_index_t otherCount = stats._otherLength.getCount() -
     stats._insertionLength.getCount() -
     stats._deletionLength.getCount() - 
     stats._inversionLength.getCount() -
     stats._duplicationLength.getCount() -
     stats._transpositionLength.getCount();
  otherCount = max(otherCount, (hal_index_t)0);
       
  os << stats._branchLength << ", " 
     << stats._genomeLength << ", " << stats._parentLength << ", "
     << stats._subs << ", "
     << stats._transitions << ", "
     << stats._transversions << ", "
     << stats._matches << ", "
     << stats._gapInsertionLength.getCount() << ", " 
     << stats._gapInsertionLength.getSum() << ", "
     << stats._gapDeletionLength.getCount() << ", " 
     << stats._gapDeletionLength.getSum() << ", "
     << stats._insertionLength.getCount() << ", " 
     << stats._insertionLength.getSum() << ", "
     << stats._deletionLength.getCount() << ", " 
     << stats._deletionLength.getSum() << ", "
     << stats._inversionLength.getCount() << ", " 
     << stats._inversionLength.getSum() << ", "
     << stats._duplicationLength.getCount() << ", " 
     << stats._duplicationLength.getSum() << ", "
     << stats._transpositionLength.getCount() << ", "
     << stats._transpositionLength.getSum() << ", "
     << otherCount;
  return os;
}
