/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <limits>
#include "halLodNode.h"
#include "halLodEdge.h"
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

  optimizeByExtension();
  optimizeByInsertion();
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
      hal_index_t maxTry = min(pos + (hal_index_t)_step / 2, 
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
  bool alreadyThere = false;
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  for (; colMapIt != colMap->end(); ++colMapIt)
  {
    const ColumnIterator::DNASet* dnaSet = colMapIt->second;
    if (dnaSet->size() > 0)
    {
      const Sequence* sequence = colMapIt->first;
      ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
      hal_index_t pos = (*dnaIt)->getArrayIndex();
      LodSegment segment(sequence, pos, false);
      SequenceMapIterator smi = _seqMap.find(sequence);
      if (smi != _seqMap.end())
      {
        SegmentSet* segmentSet = smi->second;
        SegmentIterator si = segmentSet->lower_bound(&segment);
        SegmentSet::value_compare segPLess;
        if (si != segmentSet->end() && !segPLess(&segment, *si))
        {
          alreadyThere = true;
        }
      }
      break;
    }
  }
  bool canAdd = !alreadyThere;
  hal_index_t refPos = colIt->getReferenceSequencePosition();
  if (canAdd == true && refPos != 0 && (hal_size_t)refPos != 
      colIt->getReferenceSequence()->getSequenceLength())
  {
//    canAdd = testheuristics here...
  }
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
  
  LodSegment* segment = new LodSegment(sequence, -1, false);
  _telomeres.addSegment(segment);
  segSet->insert(segment);
  segment = new LodSegment(sequence, sequence->getEndPosition() + 1, false);
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
      LodSegment* segment = new LodSegment(sequence, pos, reversed);
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
  for (BlockIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
  {
    // extend here
  }
}

void LodGraph::optimizeByInsertion()
{
  for (BlockIterator bi = _blocks.begin(); bi != _blocks.end(); ++bi)
  {
    // insert here
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
    maxSegments = max(maxSegments, (*bi)->getNumSegments());
    minSegments = min(minSegments, (*bi)->getNumSegments());

    totalLength += (*bi)->getLength();
    maxLength = max(maxLength, (*bi)->getLength());
    minLength = min(minLength, (*bi)->getLength());

    totalAdjLength += (*bi)->getTotalAdjLength();
    maxAdjLength = max(maxAdjLength, (*bi)->getTotalAdjLength());
    minAdjLength = min(minAdjLength, (*bi)->getTotalAdjLength());
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
