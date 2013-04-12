/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include <locale>
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
    if (_justSubs || i->first.second.empty() == false)
    {
      const MutationsStats& stats = i->second;       
      outStream << i->first.first << ", " << i->first.second << ", "
                << stats << endl;
    }
  }

  MutationsStats summaryStats = {0};
  hal_size_t N = 0;
  i = _branchMap.begin();
  for (; i != _branchMap.end(); ++i)
  {
    if (_justSubs || i->first.second.empty() == false)
    {
      summaryStats += i->second;
      ++N;
    }
  }
  outStream << "Total, ," << summaryStats << endl;
  if (N > 0)
  {
    summaryStats /= N;
  }
  outStream << "Average, ," << summaryStats << endl;
     
  outStream << endl;
}

void SummarizeMutations::analyzeAlignment(AlignmentConstPtr alignment,
                                          hal_size_t gapThreshold,
                                          double nThreshold, bool justSubs,
                                          const set<string>* targetSet)
{
  _gapThreshold = gapThreshold;
  _nThreshold = nThreshold;
  _justSubs = justSubs;
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
  const Genome* parent = genome->getParent();
  MutationsStats stats = {0};
  stats._genomeLength = genome->getSequenceLength();
  if (parent != NULL)
  {
    stats._parentLength = parent->getSequenceLength();
    stats._branchLength = _alignment->getBranchLength(parent->getName(),
                                                      genome->getName());
  }

  if (_justSubs == true)
  {
    substitutionAnalysis(genome, stats);
  }
  else if (parent != NULL && 
           (!_targetSet || _targetSet->find(genomeName) != _targetSet->end()))
  {
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    TopSegmentIteratorConstPtr topEnd = genome->getTopSegmentEndIterator();
    BottomSegmentIteratorConstPtr parIt = parent->getBottomSegmentIterator();

    rearrangementAnalysis(genome, stats);
  }
  
  string pname = parent != NULL ? parent->getName() : string();
  StrPair branchName(genome->getName(), pname);
  _branchMap.insert(pair<StrPair, MutationsStats>(branchName, stats));

  _alignment->closeGenome(genome);
  if (parent != NULL)
  {
    _alignment->closeGenome(parent);
  }
  vector<string> children = _alignment->getChildNames(genomeName);
  for (hal_size_t i = 0; i < children.size(); ++i)
  {
    analyzeGenomeRecursive(children[i]);
  }
}

// quickly count subsitutions without loading rearrangement machinery.
// used for benchmarks for basic file scanning... and not much else since
// the interface is still a bit wonky.
void SummarizeMutations::substitutionAnalysis(const Genome* genome, 
                                               MutationsStats& stats)
{
  assert(stats._subs == 0);
  if (genome->getNumChildren() == 0 || genome->getNumBottomSegments() == 0 ||
      (_targetSet && _targetSet->find(genome->getName()) == _targetSet->end()))
  {
    return;
  }
  const Genome* parent = genome->getParent();
  string pname = parent != NULL ? parent->getName() : string();
  StrPair branchName(genome->getName(), pname);

  BottomSegmentIteratorConstPtr bottom = genome->getBottomSegmentIterator();
  TopSegmentIteratorConstPtr top = genome->getChild(0)->getTopSegmentIterator();
  
  string gString, cString;

  hal_size_t n = genome->getNumBottomSegments();
  vector<hal_size_t> children;
  hal_size_t m = genome->getNumChildren();
  for (hal_size_t i = 0; i < m; ++i)
  {
    string cName = genome->getChild(i)->getName();
    if (!_targetSet || 
        (_targetSet && _targetSet->find(cName) != _targetSet->end()))
    {
      children.push_back(i);
    }
  }
  if (children.empty())
  {
    return;
  }

  for (hal_size_t i = 0; i < n; ++i)
  {
    bool readString = false;
    for (size_t j = 0; j < children.size(); ++j)
    {
      if (bottom->hasChild(children[j]))
      {
        if (readString == false)
        {
          bottom->getString(gString);
          readString = true;
        }
        top->toChild(bottom, children[j]);
        top->getString(cString);
        assert(gString.length() == cString.length());
        for (hal_size_t k = 0; k < gString.length(); ++k)
        {
          if (isSubstitution(gString[k], cString[k]))
          {
            ++stats._subs;
          }
        }
      }
    }
    bottom->toRight();
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
        else if (!isMissingData(child[j]) && !isMissingData(parent[j]))
        {
          ++stats._matches;
        }
      }
    }
  }
}

