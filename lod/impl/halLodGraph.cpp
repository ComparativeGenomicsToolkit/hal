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
  _genomes.clear();
  _telomeres.clear();
}

void LodGraph::build(AlignmentConstPtr alignment, const Genome* parent,
                     const vector<const Genome*>& children, 
                     hal_size_t step)
{
  erase();
  _alignment = alignment;
  _parent = parent;
  _step = step;

  assert(_parent != NULL);
  assert(_alignment->openGenome(_parent->getName()) == _parent);
  
  _genomes.insert(parent);
  for (vector<const Genome*>::const_iterator child = children.begin();
       child != children.end(); ++child)
  {
    _genomes.insert(*child);
  }
  assert(_genomes.size() == children.size() + 1);

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
  for (; seqIt != seqEnd; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    hal_size_t len = sequence->getSequenceLength();
    
    addTelomeres(sequence);
    
    for (hal_index_t pos = 0; pos < (hal_index_t)len; pos += (hal_index_t)_step)
    {
      // clamp to last position
      if (pos > 0 && pos + (hal_index_t)_step >= (hal_index_t)len)
      {
        pos =  (hal_index_t)len - 1;
      }      
      // better to move column iterator rather than getting each time?
      // -- probably not because we have to worry about dupe cache
      ColumnIteratorConstPtr colIt = 
         sequence->getColumnIterator(&_genomes, 0, pos);
      assert(colIt->getReferenceSequencePosition() == pos);
      
      // scan up to here trying to find a column we're happy to add
      hal_index_t maxTry = std::min(pos + (hal_index_t)_step / 2, 
                               (hal_index_t)len - 1);     
      hal_index_t tryPos = NULL_INDEX; 
      do 
      {
        tryPos = colIt->getReferenceSequencePosition();
        bool canAdd = canAddColumn(colIt);
        if (canAdd)
        {
          createColumn(colIt);
          break;
        }
        colIt->toRight();
      } 
      while (colIt->lastColumn() == false && tryPos < maxTry);      
    }
  }
}

bool LodGraph::canAddColumn(ColumnIteratorConstPtr colIt)
{
  // check that block has not already been added.
  hal_size_t deltaMax = 0;
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  bool breakOut = false;
  for (; colMapIt != colMap->end() && !breakOut; ++colMapIt)
  {
    const ColumnIterator::DNASet* dnaSet = colMapIt->second;
    const Sequence* sequence = colMapIt->first;
    ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
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
          deltaMax = numeric_limits<hal_size_t>::max();
          breakOut = true;
        }
        else if (!segPLess(&segment, *si))
        {
          assert(deltaMax == 0);
          breakOut = true;
        }
        else
        {
          deltaMax = std::max(deltaMax, 
                              (hal_size_t)std::abs((*si)->getLeftPos() - 
                                                   segment.getLeftPos()));
          if (si != segmentSet->begin())
          {
            --si;
            deltaMax = std::max(deltaMax, 
                                (hal_size_t)std::abs((*si)->getLeftPos() - 
                                                     segment.getLeftPos()));
          }
        }
      }
    }
  }
  bool canAdd = deltaMax > 0;
  hal_index_t refPos = colIt->getReferenceSequencePosition();
  if (canAdd == true && refPos != 0 && (hal_size_t)refPos != 
      colIt->getReferenceSequence()->getSequenceLength())
  {
    canAdd = deltaMax  >= _step / 2;

    // other heuristics here?
  }
  if (!canAdd) cout << deltaMax << "\n";
  return canAdd;
}

void LodGraph::addTelomeres(const Sequence* sequence)
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
    for (hal_size_t i = 0; i < covMapIt->second->size() && success; ++i)
    {
      if (covMapIt->second->at(i) != true)
      {
        if (success == true)
        {
          cerr << "coverage gap found at position " << i << " in"
               << " genome " << covMapIt->first->getName() << endl;
        }
        success = false;
      }
    }
    delete covMapIt->second;
  }
  return success;
}
