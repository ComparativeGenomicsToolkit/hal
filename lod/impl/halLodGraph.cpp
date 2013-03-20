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
    hal_index_t seqStart = sequence->getStartPosition();
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
      LodNode* node = new LodNode(sequence, seqStart + pos, seqStart + pos);
      _adjTable.addNode(node, colIt);
    }
  }
}
