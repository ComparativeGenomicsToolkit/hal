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

HalCons::HalCons(AlignmentConstPtr alignment, hal_size_t gapThreshold) :
  _gapThreshold(gapThreshold)
{
  analyzeAlignment(alignment);
}

HalCons::~HalCons()
{

}

void HalCons::printCsv(ostream& outStream) const
{
  outStream << "GenomeName, ParentName, BranchLength, GenomeLength," 
     " ParentLength, Subtitutions, Insertions, InsertedBases,"
     " Deletions, DeletionBases, Inversions,"
     " InvertedBases, Duplications, DuplicatedBases, Transpositions,"
     " TranspositionBases, Other, GapInsertions," 
     " GapInsertedBases, GapDeletions, GapDeletedBases" 
            << endl;

  BranchMap::const_iterator i = _branchMap.begin();
  for (; i != _branchMap.end(); ++i)
  {
    const ConsStats& stats = i->second;

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
       
    outStream << i->first.first << ", " << i->first.second << ", "
              << stats._branchLength << ", " 
              << stats._genomeLength << ", " << stats._parentLength << ", "
              << stats._subs << ", "
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
              << otherCount << ", "
              << stats._gapInsertionLength.getCount() << ", " 
              << stats._gapInsertionLength.getSum() << ", "
              << stats._gapDeletionLength.getCount() << ", " 
              << stats._gapDeletionLength.getSum() 
//              << stats._nothingLength.getCount() << ", " 
//              << stats._nothingLength.getSum()
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
  ConsStats stats = {0};

  const Genome* parent = genome->getParent();
  if (parent != NULL)
  {
    rearrangementAnalysis(genome, stats);
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    TopSegmentIteratorConstPtr topEnd = genome->getTopSegmentEndIterator();
    BottomSegmentIteratorConstPtr parIt = parent->getBottomSegmentIterator();
    string strBuf;
    string strBufParent;

    stats._genomeLength = genome->getSequenceLength();
    stats._parentLength = parent->getSequenceLength();
    stats._branchLength = _alignment->getBranchLength(parent->getName(),
                                                      genome->getName());
    
    for (; topIt != topEnd; topIt->toRight())
    {
      if (topIt->hasParent() == false || topIt->hasNextParalogy() == false)
      {
        ++stats._numInserts;
        stats._numInsertBases += topIt->getLength();
      }
      if (topIt->hasParent() == true)
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

void HalCons::rearrangementAnalysis(const Genome* genome, ConsStats& stats)
{
  const Genome* parent = genome->getParent();
  hal_index_t childIndex = parent->getChildIndex(genome);

  StrPair branchName(genome->getName(), parent->getName());

  // do the gapped deletions by scanning the parent
  GappedBottomSegmentIteratorConstPtr gappedBottom = 
     parent->getGappedBottomSegmentIterator(0, childIndex, _gapThreshold);

  while (gappedBottom->getRightArrayIndex() < parent->getNumBottomSegments())
  {
    hal_size_t numGaps = gappedBottom->getNumGaps();
    if (numGaps > 0)
    {
      stats._gapDeletionLength.add(gappedBottom->getNumGapBases(), numGaps);
    }
    gappedBottom->toRight();
  }

  GappedTopSegmentIteratorConstPtr gappedTop = 
     genome->getGappedTopSegmentIterator(0, _gapThreshold);

  while (gappedTop->getRightArrayIndex() < genome->getNumTopSegments())
  {
    hal_size_t numGaps = gappedTop->getNumGaps();
    if (numGaps > 0)
    {
      stats._gapInsertionLength.add(gappedTop->getNumGapBases(), numGaps);
    }
    gappedTop->toRight();
  }

  RearrangementPtr r = genome->getRearrangement();
  r->setGapLengthThreshold(_gapThreshold);
  do {
    switch (r->getID())
    {
    case Rearrangement::Inversion:
      stats._inversionLength.add(r->getLength());
      break;
    case Rearrangement::Deletion:
      stats._deletionLength.add(r->getLength());
      break;
    case Rearrangement::Transposition:
      stats._transpositionLength.add(r->getLength());
      break;
    case Rearrangement::Insertion:
      stats._insertionLength.add(r->getLength());
      break;
    case Rearrangement::Duplication:
      stats._duplicationLength.add(r->getLength());
      break;
    case Rearrangement::Nothing:
      stats._nothingLength.add(r->getLength());
      break;
    default:
      stats._otherLength.add(r->getLength());
      break;
    }
  } 
  while (r->identifyNext() == true);
}
