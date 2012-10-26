/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halSummarizeMutations.h"

using namespace std;
using namespace hal;

SummarizeMutations::SummarizeMutations()
{

}

SummarizeMutations::~SummarizeMutations()
{

}

void SummarizeMutations::printCsv(ostream& outStream) const
{
  outStream << "GenomeName, ParentName, ";
  MutationsStats::printHeader(outStream);
  outStream << endl;

  BranchMap::const_iterator i = _branchMap.begin();
  for (; i != _branchMap.end(); ++i)
  {
    const MutationsStats& stats = i->second;       
    outStream << i->first.first << ", " << i->first.second << ", "
              << stats << endl;
  }

  outStream << endl;
}

void SummarizeMutations::analyzeAlignment(AlignmentConstPtr alignment,
                                          hal_size_t gapThreshold,
                                          double nThreshold,
                                          const set<string>* targetSet)
{
  _gapThreshold = gapThreshold;
  _nThreshold = nThreshold;
  _targetSet = targetSet;
  _branchMap.clear();
  _alignment = alignment;
  
  if (_alignment->getNumGenomes() > 0)
  {
    string root = _alignment->getRootName();
    analyzeGenomeRecursive(root);
  }

  _alignment = AlignmentConstPtr();
}

void SummarizeMutations::analyzeGenomeRecursive(const string& genomeName)
{
  const Genome* genome = _alignment->openGenome(genomeName);
  assert(genome != NULL);
  MutationsStats stats = {0};

  const Genome* parent = genome->getParent();
  if (parent != NULL && 
      (!_targetSet || _targetSet->find(genomeName) != _targetSet->end()))
  {
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    TopSegmentIteratorConstPtr topEnd = genome->getTopSegmentEndIterator();
    BottomSegmentIteratorConstPtr parIt = parent->getBottomSegmentIterator();
    string strBuf;
    string strBufParent;

    stats._genomeLength = genome->getSequenceLength();
    stats._parentLength = parent->getSequenceLength();
    stats._branchLength = _alignment->getBranchLength(parent->getName(),
                                                      genome->getName());
   
    rearrangementAnalysis(genome, stats);
     
    StrPair branchName(genome->getName(), parent->getName());
    _branchMap.insert(pair<StrPair, MutationsStats>(branchName, stats));

    _alignment->closeGenome(parent);
  }
  _alignment->closeGenome(genome);
  vector<string> children = _alignment->getChildNames(genomeName);
  for (hal_size_t i = 0; i < children.size(); ++i)
  {
    analyzeGenomeRecursive(children[i]);
  }
}

void SummarizeMutations::rearrangementAnalysis(const Genome* genome, 
                                               MutationsStats& stats)
{
  const Genome* parent = genome->getParent();
  hal_index_t childIndex = parent->getChildIndex(genome);

  StrPair branchName(genome->getName(), parent->getName());

  // do the gapped deletions by scanning the parent
  GappedBottomSegmentIteratorConstPtr gappedBottom = 
     parent->getGappedBottomSegmentIterator(0, childIndex, _gapThreshold);

  while (gappedBottom->getRightArrayIndex() < 
         (hal_index_t)parent->getNumBottomSegments())
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

  RearrangementPtr r = genome->getRearrangement(0, _gapThreshold, _nThreshold);
  do {    
    // get the number of gaps from the current range of the rearrangement
    // (this should cover the entire genome)
    gappedTop->setLeft(r->getLeftBreakpoint());
    subsAndGapInserts(gappedTop, stats);
        
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

void SummarizeMutations::subsAndGapInserts(
  GappedTopSegmentIteratorConstPtr gappedTop, MutationsStats& stats)
{
  assert(gappedTop->getReversed() == false);
  hal_size_t numGaps = gappedTop->getNumGaps();
  if (numGaps > 0)
  {
    stats._gapInsertionLength.add(gappedTop->getNumGapBases(), numGaps);
  }

  string parent, child;
  TopSegmentIteratorConstPtr l = gappedTop->getLeft();
  TopSegmentIteratorConstPtr r = gappedTop->getRight();
  BottomSegmentIteratorConstPtr p = 
     l->getTopSegment()->getGenome()->getParent()->getBottomSegmentIterator();

  for (TopSegmentIteratorConstPtr i = l->copy(); 
       i->getTopSegment()->getArrayIndex() <= 
          r->getTopSegment()->getArrayIndex();
       i->toRight())
  {
    if (i->hasParent())
    {
      p->toParent(i);
      i->getString(child);
      p->getString(parent);
      assert(child.length() == parent.length());
      for (size_t j = 0; j < child.length(); ++j)
      {
        if (isTransition(child[j], parent[j]))
        {
          ++stats._transitions;
          ++stats._subs;
        }
        else if (isTransversion(child[j], parent[j]))
        {
          ++stats._transversions;
          ++stats._subs;
        }
        else if (isSubstitution(child[j], parent[j]))
        {
          ++stats._subs;
        }
      }
    }
  }
}
