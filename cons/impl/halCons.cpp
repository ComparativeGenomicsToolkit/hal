/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halCons.h"

using namespace std;
using namespace hal;

HalCons::HalCons()
{

}

HalCons::HalCons(AlignmentConstPtr alignment)
{
  analyzeAlignment(alignment);
}

HalCons::~HalCons()
{

}

void HalCons::printCsv(ostream& outStream) const
{
  outStream << "GenomeName, ParentName, BranchLength, GenomeLength," 
     " ParentLength, Subtitutions, Insertions, InsertedBases, Inversions,"
     " InvertedBases, Duplications, DuplicatedBases, GapInsertions," 
     " GapInsertedBases, GapDeletions, GapDeletedBases" << endl;

  BranchMap::const_iterator i = _branchMap.begin();
  for (; i != _branchMap.end(); ++i)
  {
    const ConsStats& stats = i->second;
    outStream << i->first.first << ", " << i->first.second << ", "
              << stats._branchLength << ", " 
              << stats._genomeLength << ", " << stats._parentLength << ", "
              << stats._subs << ", "
              << stats._numInserts << ", " << stats._numInsertBases << ", "
              << stats._numInverts << ", " << stats._numInvertBases << ", "
              << stats._numDups << ", " << stats._numDupBases << ", "
              << stats._gapInserts << ", " << stats._gapInsertBases << ", "
              << stats._gapDeletes << ", " << stats._gapDeleteBases
              <<endl;
  }

  outStream << endl;
}

void HalCons::analyzeAlignment(AlignmentConstPtr alignment)
{
  _branchMap.clear();
  _alignment = alignment;
  
  if (_alignment->getNumGenomes() > 0)
  {
    string root = _alignment->getRootName();
    analyzeGenomeRecursive(root);
  }

  _alignment = AlignmentConstPtr();
}

void HalCons::analyzeGenomeRecursive(const string& genomeName)
{
  const Genome* genome = _alignment->openGenome(genomeName);
  assert(genome != NULL);

  const Genome* parent = genome->getParent();
  if (parent != NULL)
  {
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    TopSegmentIteratorConstPtr topEnd = genome->getTopSegmentEndIterator();
    BottomSegmentIteratorConstPtr parIt = parent->getBottomSegmentIterator();
    string strBuf;
    string strBufParent;

    ConsStats stats = {0};
    stats._genomeLength = genome->getSequenceLength();
    stats._parentLength = parent->getSequenceLength();
    stats._branchLength = _alignment->getBranchLength(parent->getName(),
                                                      genome->getName());
    
    for (; topIt != topEnd; topIt->toRight())
    {
      if (topIt->getTopSegment()->isGapInsertion() == true)
      {
        ++stats._gapInserts;
        stats._gapInsertBases += topIt->getLength();
      }
      else if (topIt->hasParent() == false || topIt->hasNextParalogy() == false)
      {
        ++stats._numInserts;
        stats._numInsertBases += topIt->getLength();
      }
      else
      {
        parIt->toParent(topIt);
        if (parIt->getReversed() == true)
        {
          ++stats._numInverts;
          stats._numInvertBases += topIt->getLength();
        }
        parIt->getString(strBufParent);
        topIt->getString(strBuf);
        stats._subs += hammingDistance(strBuf, strBufParent);
        
        if (topIt->hasNextParalogy() == true)
        {
          ++stats._numDups;
          stats._numDupBases += topIt->getLength();
        }
      }
    }
    //TODO: DELETIONS!!!
    
    StrPair branchName(genome->getName(), parent->getName());
    _branchMap.insert(pair<StrPair, ConsStats>(branchName, stats));

    _alignment->closeGenome(parent);
  }
  _alignment->closeGenome(genome);
  vector<string> children = _alignment->getChildNames(genomeName);
  for (hal_size_t i = 0; i < children.size(); ++i)
  {
    analyzeGenomeRecursive(children[i]);
  }
}
