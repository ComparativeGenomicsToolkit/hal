/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include "halLodNode.h"
#include "halLodEdge.h"
#include "halLodGraph.h"

using namespace std;
using namespace hal;

LodGraph::LodGraph()
{

}

LodGraph::~LodGraph()
{
  erase();
}

void LodGraph::erase()
{
  for (GenomeNodesIterator i = _genomeNodes.begin(); i != _genomeNodes.end();
       ++i)
  {
    NodeList* nodeList = i->second;
    for (NodeIterator j = nodeList->begin(); j != nodeList->end(); ++j)
    {
      delete *j;
    }
    delete nodeList;
    _genomeNodes.erase(i);
  }
  _children.clear();
  _parent = NULL;
  _alignment = AlignmentConstPtr();
}

void LodGraph::build(AlignmentConstPtr alignment, const Genome* parent,
                     const vector<const Genome*> children, 
                     hal_size_t step)
{
  erase();
  _alignment = alignment;
  _parent = parent;
  _children = children;
  _step = step;

  assert(_parent != NULL);
  assert(_alignment->openGenome(_parent->getName()) == _parent);
  
  _genomeNodes.insert(pair<const Genome*, NodeList*>(_parent, new NodeList()));
  for (size_t i = 0; i < _children.size(); ++i)
  {
    _genomeNodes.insert(pair<const Genome*, NodeList*>(_children[i], 
                                                       new NodeList()));
  }

  for (GenomeNodesIterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    scanGenome(gni->first, gni->second);
  }

  _adjTable.writeAdjacenciesIntoNodes();
  _adjTable.clear();
}


void LodGraph::scanGenome(const Genome* genome, NodeList* nodeList)
{
  // should be done at outside scope.
  set<const Genome*> tgtSet;
  for (GenomeNodesIterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    tgtSet.insert(gni->first);
  }

  SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
  SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
  for (; seqIt != seqEnd; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    hal_size_t len = sequence->getSequenceLength();

    for (hal_index_t pos = 0; pos < (hal_index_t)len; pos += (hal_index_t)_step)
    {
      // clamp to last position
      if (pos > 0 && pos + (hal_index_t)_step >= (hal_index_t)len)
      {
        pos =  (hal_index_t)len - 1;
      }      
      // better to move column iterator rather than getting each time?
      ColumnIteratorConstPtr colIt = sequence->getColumnIterator(&tgtSet, 0,
                                                                 pos);
      // convert pos to genome coordinate
      createColumn(colIt, nodeList);
    }
  }
}

void LodGraph::createColumn(ColumnIteratorConstPtr colIt, NodeList* nodeList)
{
  const Genome* refGenome = colIt->getReferenceGenome();
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  for (; colMapIt != colMap->end(); ++colMapIt)
  {
    // add a new node for every palagous base in the reference genome
    if (colMapIt->first->getGenome() == refGenome)
    {
      const ColumnIterator::DNASet* dnaSet = colMapIt->second;
      const Sequence* sequence = colMapIt->first;
      ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
      for (; dnaIt != dnaSet->end(); ++dnaIt)
      {
        LodNode* node = new LodNode(sequence, (*dnaIt)->getArrayIndex(), 
                                    (*dnaIt)->getArrayIndex());
        bool added = _adjTable.addNode(node, colIt);
        if (added == true)
        {
          nodeList->push_back(node);
        }
        else
        {
          delete node;
        }
      }      
    }
  }
}
