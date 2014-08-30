/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "halLodGraph.h"

using namespace std;
using namespace hal;

LodGraph::LodGraph() : _extendFraction(1.0)
{

}

LodGraph::~LodGraph()
{
  erase();
}

void LodGraph::erase()
{
  for (SequenceMapIterator smi = _seqMap.begin(); smi != _seqMap.end(); ++smi)
  {
    delete smi->second;
  }
  _seqMap.clear();
  for (BlockIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
  {
    delete *bi;
  }
  _blocks.clear();
  _parent = NULL;
  _grandParent = NULL;
  _genomes.clear();
  _telomeres.clear();
}

void LodGraph::build(AlignmentConstPtr alignment, const Genome* parent,
                     const vector<const Genome*>& children, 
                     const Genome* grandParent,
                     hal_size_t step, bool allSequences, double probeFrac,
                     double minSeqFrac)
{
  erase();
  _alignment = alignment;
  _parent = parent;
  _grandParent = grandParent;
  _step = step;
  _allSequences = allSequences;
  _probeFrac = probeFrac;
  _minSeqLen = minSeqFrac * _step;

  assert(_parent != NULL);
  assert(_alignment->openGenome(_parent->getName()) == _parent);
  
  _genomes.insert(parent);
  for (vector<const Genome*>::const_iterator child = children.begin();
       child != children.end(); ++child)
  {
    _genomes.insert(*child);
  }
  assert(_genomes.size() == children.size() + 1);
  if (_grandParent != NULL)
  {
    _genomes.insert(_grandParent);
  }

  for (set<const Genome*>::iterator gi = _genomes.begin(); 
       gi != _genomes.end(); ++gi)
  {
    scanGenome(*gi);
  }
  
  computeAdjacencies();
  printDimensions(cout);
  optimizeByExtension();
  printDimensions(cout);
  optimizeByMerging();
  printDimensions(cout);
  optimizeByInsertion();
  printDimensions(cout);
  assert(checkCoverage() == true);
}


void LodGraph::scanGenome(const Genome* genome)
{
  SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
  SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
  hal_index_t lastSampledPos = 0;
  hal_index_t halfStep = std::max((hal_index_t)1, (hal_index_t)_step / 2);
  for (; seqIt != seqEnd; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    hal_size_t len = sequence->getSequenceLength();
    hal_index_t seqEnd = sequence->getStartPosition() + (hal_index_t)len;

    addTelomeres(sequence);
    if (_allSequences == true ||
        (sequence->getSequenceLength() > _minSeqLen &&
         seqEnd - lastSampledPos > (hal_index_t)_step))
    {
      for (hal_index_t pos = 0; pos < (hal_index_t)len; 
           pos += (hal_index_t)_step)
      {
        // clamp to last position
        if (pos > 0 && pos + (hal_index_t)_step >= (hal_index_t)len)
        {
          pos =  (hal_index_t)len - 1;
        }   

        // scan range trying to find genome to add
        hal_index_t minTry = std::max((hal_index_t)0, pos - halfStep);
        hal_index_t maxTry = std::min(pos + halfStep, (hal_index_t)len - 1);
        double redProbFac = _probeFrac * (
          (double)(maxTry - minTry) / (double)sequence->getSequenceLength());
        hal_index_t numProbe = 
           (hal_index_t)std::max(1., (double)(maxTry - minTry) * redProbFac);
        hal_index_t npMinus1 = numProbe < 2 ? numProbe : numProbe - 1;
        hal_index_t probeStep = std::max((hal_index_t)1,
                                         (maxTry - minTry) / (npMinus1));
        hal_index_t bestPos = NULL_INDEX;
        hal_size_t maxNumGenomes = 1;
        hal_size_t maxDelta = 0;     
        hal_size_t maxMinSeqLen = 0;
        hal_index_t tryPos = numProbe == 1 ? pos : minTry;
        ColumnIteratorConstPtr colIt = 
           sequence->getColumnIterator(&_genomes, 0, tryPos);
        do 
        {
          if (colIt->getReferenceSequencePosition() != tryPos)
          {
            colIt->toSite(sequence->getStartPosition() + tryPos, 
                          sequence->getEndPosition(), true);
          }
          assert(colIt->getReferenceSequence() == sequence);
          assert(colIt->getReferenceSequencePosition() == tryPos);
          hal_size_t delta;
          hal_size_t numGenomes;
          hal_size_t minSeqLen;
          evaluateColumn(colIt, delta, numGenomes, minSeqLen);
          if (bestColumn(probeStep, delta, numGenomes, minSeqLen,
                         maxDelta, maxNumGenomes, maxMinSeqLen))
          {
            bestPos = tryPos;
            maxDelta = delta;
            maxNumGenomes = numGenomes;
            maxMinSeqLen = minSeqLen;
          }
          tryPos += probeStep;
        } 
        while (colIt->lastColumn() == false && tryPos < maxTry);

        if (bestPos != NULL_INDEX)
        {
          if (colIt->getReferenceSequencePosition() != bestPos)
          {
            colIt->toSite(sequence->getStartPosition() + bestPos, 
                          sequence->getEndPosition(), true);
          }
          assert(colIt->getReferenceSequence() == sequence);
          assert(colIt->getReferenceSequencePosition() == bestPos);
          createColumn(colIt);
          lastSampledPos = sequence->getStartPosition() + bestPos;
        }
      }
    }
  }
}

void LodGraph::evaluateColumn(ColumnIteratorConstPtr colIt,
                              hal_size_t& outDeltaMax,
                              hal_size_t& outNumGenomes,
                              hal_size_t& outMinSeqLen)
{
  outDeltaMax = 0;
  outMinSeqLen = numeric_limits<hal_size_t>::max();
  set<const Genome*> genomeSet;
  // check that block has not already been added.
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  bool breakOut = false;
  for (; colMapIt != colMap->end() && !breakOut; ++colMapIt)
  {
    const ColumnIterator::DNASet* dnaSet = colMapIt->second;
    const Sequence* sequence = colMapIt->first;
    if (sequence->getSequenceLength() <= _minSeqLen)
    {
      // we never want to align two leaves through a disappeared
      // contig in parent
      if (sequence->getGenome() == _parent)
      {
        outDeltaMax = 0;
        outNumGenomes = 0;
        outMinSeqLen = 0;
        break;
      }
    }
    else
    {
      outMinSeqLen = std::min(outMinSeqLen, sequence->getSequenceLength());
      ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
      if (!dnaSet->empty())
      {
        genomeSet.insert(sequence->getGenome());
      }
      for (; dnaIt != dnaSet->end() && !breakOut; ++dnaIt)
      {
        hal_index_t pos = (*dnaIt)->getArrayIndex();
        LodSegment segment(NULL, sequence, pos, false);
        SequenceMapIterator smi = _seqMap.find(sequence);
        if (smi != _seqMap.end())
        {
          SegmentSet* segmentSet = smi->second;
          SegmentIterator si = segmentSet->lower_bound(&segment);
          SegmentSet::value_compare segPLess = segmentSet->key_comp();
          if (si == segmentSet->end())
          {
            outDeltaMax = numeric_limits<hal_size_t>::max();
            breakOut = true;
          }
          else if (!segPLess(&segment, *si))
          {
            assert(outDeltaMax == 0);
            breakOut = true;
          }
          else
          {
            hal_size_t delta = 
               std::min(_step, (hal_size_t)std::abs((*si)->getLeftPos() - 
                                                    segment.getLeftPos()));
            if (si != segmentSet->begin())
            {
              --si;
              delta += (hal_size_t)std::abs((*si)->getLeftPos() - 
                                            segment.getLeftPos());
            }
            else
            {
              delta *= 2;
            }
            outDeltaMax = std::max(outDeltaMax, delta);
          }
        }
      }
    }
  }
  outNumGenomes = genomeSet.size();
}

bool LodGraph::bestColumn(hal_size_t probeStep, hal_size_t delta, 
                          hal_size_t numGenomes, hal_size_t minSeqLen,
                          hal_size_t maxDelta, hal_size_t maxNumGenomes,
                          hal_size_t maxMinSeqLen)
{
  if (delta <= (hal_size_t)probeStep)
  {
    return false;
  }
  if (numGenomes > maxNumGenomes)
  {
    return true;
  }
  else if (numGenomes == maxNumGenomes)
  {
    if (minSeqLen > maxMinSeqLen)
    {
      return true;
    }
    else if (minSeqLen == maxMinSeqLen)
    {
      return numGenomes > 1 && delta > maxDelta;
    }
  }
  return false;
}

void LodGraph::addTelomeres(const Sequence* sequence)
{
  SequenceMapIterator smi = _seqMap.find(sequence);
  SegmentSet* segSet = NULL;
  if (smi == _seqMap.end())
  {
    segSet = new SegmentSet();
    _seqMap.insert(pair<const Sequence*, SegmentSet*>(sequence, segSet));
  }
  else
  {
    segSet = smi->second;
  }

  LodSegment* segment = new LodSegment(&_telomeres, sequence, 
                                       sequence->getStartPosition() - 1,
                                       false);
  _telomeres.addSegment(segment);
  segSet->insert(segment);  
  segment = new LodSegment(&_telomeres, sequence, 
                           sequence->getEndPosition() + 1, false);
  _telomeres.addSegment(segment);
  segSet->insert(segment);
}

void LodGraph::createColumn(ColumnIteratorConstPtr colIt)
{
  LodBlock* block = new LodBlock();
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  for (; colMapIt != colMap->end(); ++colMapIt)
  {
    const Sequence* sequence = colMapIt->first;
    if (sequence->getSequenceLength() > _minSeqLen)
    {
      SequenceMapIterator smi = _seqMap.find(sequence);
      SegmentSet* segSet = NULL;
      if (smi == _seqMap.end())
      {
        segSet = new SegmentSet();;
        _seqMap.insert(pair<const Sequence*, SegmentSet*>(sequence, segSet));
      }
      else
      {
        segSet = smi->second;
      }
    
      const ColumnIterator::DNASet* dnaSet = colMapIt->second;
      for (ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
           dnaIt != dnaSet->end(); ++ dnaIt)
      {
        hal_index_t pos = (*dnaIt)->getArrayIndex();
        bool reversed = (*dnaIt)->getReversed();
        LodSegment* segment = new LodSegment(block, sequence, pos, reversed);
        block->addSegment(segment);
        assert(segSet->find(segment) == segSet->end());
        segSet->insert(segment);
      }
    }
  }
  assert(block->getNumSegments() > 0);
  _blocks.push_back(block);
}

void LodGraph::computeAdjacencies()
{
  for (SequenceMapIterator smi = _seqMap.begin(); smi != _seqMap.end(); ++smi)
  {
    SegmentIterator si = smi->second->begin();
    SegmentIterator siNext;
    for (; si != smi->second->end(); ++si)
    {
      siNext = si;
      ++siNext;
      if (siNext != smi->second->end())
      {
        assert((*si)->overlaps(**siNext) == false);
        (*si)->addEdgeFromRightToLeft(*siNext);
      }
    }
  }
}

void LodGraph::optimizeByExtension()
{
  // Put bigger blocks first because (we hope) they tend to represent
  // more information and we'd rather extend them than, say, a block
  // that represents a little insertion. 
  std::sort(_blocks.begin(), _blocks.end(), LodBlockPBigger());
  
  // Pass 1: Extend all blocks with at least 2 sequences by half
  for (BlockIterator bi = _blocks.begin(); bi != _blocks.end() && 
          (*bi)->getNumSegments() > 1; ++bi)
  {
    (*bi)->extend(0.5);
  }

  // Pass 2: Greedily extend each block to the max
  for (BlockIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
  {
    (*bi)->extend();
  }

}

void LodGraph::optimizeByMerging()
{
  BlockList mergeList(_blocks);
  for (BlockIterator bi = mergeList.begin(); bi != mergeList.end(); ++bi)
  {
    LodBlock* adjBlock = (*bi)->getHeadMergePartner();
    if (adjBlock != NULL)
    {
      for (hal_size_t segIdx = 0; segIdx < adjBlock->getNumSegments(); ++segIdx)
      {
        // blah - need to clean interface but this is harmless for now
        LodSegment* seg = const_cast<LodSegment*>(
          adjBlock->getSegment(segIdx));
        SequenceMapIterator mapIt = _seqMap.find(seg->getSequence());
        assert(mapIt != _seqMap.end());
        SegmentIterator setIt = mapIt->second->find(seg);
        assert(setIt != mapIt->second->end());
        assert(*setIt == seg);
        mapIt->second->erase(setIt);
      }
      (*bi)->mergeHead(adjBlock);
    }
  }
  _blocks.clear();
  for (BlockIterator bi = mergeList.begin(); bi != mergeList.end(); ++bi)
  {
    if ((*bi)->getNumSegments() > 0)
    {
      _blocks.push_back(*bi);
    }
  }
}

void LodGraph::optimizeByInsertion()
{
  vector<LodBlock*> newBlocks;
  BlockIterator startPoint = _blocks.begin();
  // seems more convoluted than necessary but I had problems with ?iterators?
  // doing it more simply. 
  while (startPoint != _blocks.end())
  {
    newBlocks.clear();
    for (BlockIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
    {
      (*bi)->insertNeighbours(newBlocks);
    }
    startPoint = _blocks.end();
    --startPoint;
    _blocks.insert(_blocks.end(), newBlocks.begin(), newBlocks.end());

    // need to get the new segments into the sorted structure too!
    for (BlockIterator bi = newBlocks.begin(); bi != newBlocks.end(); ++bi)
    {
      for (hal_size_t i = 0; i < (*bi)->getNumSegments(); ++i)
      {
        // blah - need to clean interface but this is harmless for now
        LodSegment* seg = const_cast<LodSegment*>((*bi)->getSegment(i));
        _seqMap.find(seg->getSequence())->second->insert(seg);
      }
    }
    ++startPoint;
  }
}

void LodGraph::printDimensions(ostream& os) const
{
  hal_size_t totalSegments = 0;
  hal_size_t maxSegments = 0;
  hal_size_t minSegments = numeric_limits<hal_size_t>::max();

  hal_size_t totalLength = 0;
  hal_size_t maxLength = 0;
  hal_size_t minLength = numeric_limits<hal_size_t>::max();

  hal_size_t totalAdjLength = 0;
  hal_size_t maxAdjLength = 0;
  hal_size_t minAdjLength = numeric_limits<hal_size_t>::max();

  for (BlockConstIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
  {
    totalSegments += (*bi)->getNumSegments();
    maxSegments = std::max(maxSegments, (*bi)->getNumSegments());
    minSegments = std::min(minSegments, (*bi)->getNumSegments());

    totalLength += (*bi)->getLength();
    maxLength = std::max(maxLength, (*bi)->getLength());
    minLength = std::min(minLength, (*bi)->getLength());

    totalAdjLength += (*bi)->getTotalAdjLength();
    maxAdjLength = std::max(maxAdjLength, (*bi)->getTotalAdjLength());
    minAdjLength = std::min(minAdjLength, (*bi)->getTotalAdjLength());
  }

  os << "Graph: numBlocks=" << _blocks.size()
     << " numSegs=" << totalSegments << " minSegs=" << minSegments
     << " maxSegs=" << maxSegments
     << " totLen=" << totalLength << " minLen=" << minLength
     << " maxLen=" << maxLength 
     << " totAdjLen=" << totalAdjLength << " minAdjLen=" << minAdjLength
     << " maxAdjLen=" << maxAdjLength 
     << endl;
}

bool LodGraph::checkCoverage() const
{
  map<const Genome*, vector<bool>* > covMap;
  pair<const Genome*, vector<bool>* > covPair;
  map<const Genome*, vector<bool>* >::iterator covMapIt;
  set<const Genome*>::const_iterator genomeIt;
  bool success = true;

  // create new bit vector for each genome, set all zero
  for (genomeIt = _genomes.begin(); genomeIt != _genomes.end(); ++genomeIt)
  {
    covPair.first = *genomeIt;
    covPair.second = new vector<bool>(covPair.first->getSequenceLength(), 
                                      false);
    covMap.insert(covPair);
  }

  // iterate through every segment in every block.
  for (BlockConstIterator bi = _blocks.begin(); bi != _blocks.end() && success;
       ++bi)
  {
    for (hal_size_t segIdx = 0; segIdx < (*bi)->getNumSegments() && success; 
         ++ segIdx)
    {
      const LodSegment* segment = (*bi)->getSegment(segIdx);
      const Genome* genome = segment->getSequence()->getGenome();
      covMapIt = covMap.find(genome);
      assert(covMapIt != covMap.end());
      assert(segment->getLeftPos() <= segment->getRightPos());
      for (hal_index_t pos = segment->getLeftPos();
           pos <= segment->getRightPos(); ++pos)
      {
        assert(pos < (hal_index_t)genome->getSequenceLength());
        if (covMapIt->second->at(pos) == true)
        {
          cerr << "dupcliate coverage found at position " << pos << " in"
               << " genome " << genome->getName() << " offending block:\n"
               << **bi << endl;
          success = false;
        }
        covMapIt->second->at(pos) = true;
      }
    }
  }

  // check the coverage and clean up
  for (covMapIt = covMap.begin(); covMapIt != covMap.end(); ++covMapIt)
  {
    SequenceIteratorConstPtr seqIt = covMapIt->first->getSequenceIterator();
    SequenceIteratorConstPtr seqEnd = covMapIt->first->getSequenceEndIterator();
    for (; seqIt != seqEnd && success; seqIt->toNext())
    {
      const Sequence* sequence = seqIt->getSequence();
      if (sequence->getSequenceLength() > 0)
      {
        hal_size_t numCovered = 0;
        for (hal_index_t i = sequence->getStartPosition(); 
             i <= sequence->getEndPosition(); ++i)
        {
          if (covMapIt->second->at(i) != true)
          {
            ++numCovered;
          }
        }
        if (numCovered != 0 && numCovered != sequence->getSequenceLength())
        {
          cerr << numCovered << " out of " << sequence->getSequenceLength()
               << " bases covered for sequence " << sequence->getFullName()
               << endl;
          success = false;
        }
        else if (numCovered == 0 && _allSequences == true)
        {
          cerr << sequence->getFullName() << " uncovered even though "
               << "allSequences set to true" << endl;
          success = false;
        }
      }
    }
    delete covMapIt->second;
  }
  return success;
}
