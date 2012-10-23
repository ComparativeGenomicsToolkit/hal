/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <deque>
#include <cassert>
#include "halMutations.h"

using namespace std;
using namespace hal;

Mutations::Mutations(AlignmentConstPtr alignment, hal_size_t gapThreshold) :
  _gapThreshold(gapThreshold)
{
  analyzeAlignment(alignment);
}

Mutations::~Mutations()
{

}

void Mutations::printCsv(ostream& outStream) const
{
  outStream << "GenomeName, ParentName, BranchLength, GenomeLength," 
     " ParentLength, Subtitutions, Transitions, Transversions,"
     " GapInsertions, GapInsertedBases, GapDeletions, GapDeletedBases," 
     " Insertions, InsertionBases, Deletions, DeletionBases, Inversions,"
     " InvertedBases, Duplications, DuplicatedBases, Transpositions,"
     " TranspositionBases, Other" 
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
              << stats._transitions << ", "
              << stats._transversions << ", "
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
              << otherCount << ", "
              <<endl;
  }

  outStream << endl;
}

void Mutations::analyzeAlignment(AlignmentConstPtr alignment)
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

void Mutations::analyzeGenomeRecursive(const string& genomeName)
{
  const Genome* genome = _alignment->openGenome(genomeName);
  assert(genome != NULL);
  ConsStats stats = {0};

  const Genome* parent = genome->getParent();
  if (parent != NULL)
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

void Mutations::rearrangementAnalysis(const Genome* genome, ConsStats& stats)
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

  RearrangementPtr r = genome->getRearrangement();
  r->setGapLengthThreshold(_gapThreshold);
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

void Mutations::subsAndGapInserts(GappedTopSegmentIteratorConstPtr gappedTop, 
                                  ConsStats& stats)
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
       i->getTopSegment()->getArrayIndex() <= r->getTopSegment()->getArrayIndex();
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
