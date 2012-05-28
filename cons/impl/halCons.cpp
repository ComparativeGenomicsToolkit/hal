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
     " InvertedBases" << endl;

  BranchMap::const_iterator i = _branchMap.begin();
  for (; i != _branchMap.end(); ++i)
  {
    const ConsStats& stats = i->second;
    outStream << i->first.first << ", " << i->first.second << ", "
              << stats._branchLength << ", " 
              << stats._genomeLength << ", " << stats._parentLength << ", "
              << stats._subs << ", "
              << stats._numInserts << ", " << stats._numInsertBases << ", "
              << stats._numInverts << ", " << stats._numInvertBases <<endl;
  }

  outStream << endl;
}

void HalCons::analyzeAlignment(AlignmentConstPtr alignment)
{
  _branchMap.clear();
  _alignment = alignment;
  
  if (_alignment->getNumGenomes() > 0)
  {
    const Genome* root = _alignment->openGenome(_alignment->getRootName());
    analyzeGenomeRecursive(root);
  }

  _alignment = AlignmentConstPtr();
}

void HalCons::analyzeGenomeRecursive(const Genome* genome)
{
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
      if (topIt->hasParent() == false)
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
      }
    }
    //TODO: DELETIONS!!!
    
    StrPair branchName(genome->getName(), parent->getName());
    _branchMap.insert(pair<StrPair, ConsStats>(branchName, stats));
  }

  for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
  {
    const Genome* child = genome->getChild(i);
    analyzeGenomeRecursive(child);
  }
}
