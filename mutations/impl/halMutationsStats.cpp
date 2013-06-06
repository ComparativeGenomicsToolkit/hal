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
     " ParentLength, Subtitutions, Transitions, Transversions, Matches,"
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

MutationsStats& hal::operator+=(MutationsStats& ms, const MutationsStats& other)
{
   ms._genomeLength += other._genomeLength;
   ms._parentLength += other._parentLength;
   ms._branchLength += other._branchLength;

   ms._subs += other._subs;
   ms._transitions += other._transitions;
   ms._transversions += other._transversions;
   ms._matches += other._matches;

   ms._nothingLength += other._nothingLength;
   ms._inversionLength += other._inversionLength;
   ms._insertionLength += other._insertionLength;
   ms._deletionLength += other._deletionLength;
   ms._transpositionLength += other._transpositionLength;
   ms._duplicationLength += other._duplicationLength;
   ms._otherLength += other._otherLength;
   ms._gapInsertionLength += other._gapInsertionLength;
   ms._gapDeletionLength += other._gapDeletionLength;

   return ms;
}

MutationsStats& hal::operator/=(MutationsStats& ms, hal_size_t N)
{
   ms._genomeLength /= N;
   ms._parentLength /= N;
   ms._branchLength /= N;

   ms._subs /= N;
   ms._transitions /= N;
   ms._transversions /= N;
   ms._matches /= N;

   ms._nothingLength /= N;
   ms._inversionLength /= N;
   ms._insertionLength /= N;
   ms._deletionLength /= N;
   ms._transpositionLength /= N;
   ms._duplicationLength /= N;
   ms._otherLength /= N;
   ms._gapInsertionLength /= N;
   ms._gapDeletionLength /= N;

   return ms;
}
