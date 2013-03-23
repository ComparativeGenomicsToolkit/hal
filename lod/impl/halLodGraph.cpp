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

  optimizeByExtension();
  optimizeByInsertion();
  
  for (GenomeNodesIterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    gni->second->sort(LodNodePLess());
  }
  
  assert(checkCoverage() == true);
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
    hal_index_t seqStart = sequence->getStartPosition();
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
      assert(colIt->getReferenceSequencePosition() == pos);
      
      // scan up to here trying to find a column we're happy to add
      hal_index_t maxTry = min(pos + (hal_index_t)_step / 2, 
                               (hal_index_t)len - 1);     
      hal_index_t tryPos = NULL_INDEX; 
      do 
      {
        tryPos = colIt->getReferenceSequencePosition();
        bool canAdd = _adjTable.canAddColumn(colIt, _step);
        if (canAdd)
        {
          createColumn(colIt, nodeList);
          break;
        }
        else if (tryPos == 0 || tryPos == (hal_index_t)len - 1)
        {
          // no choice but to add so we dump in without homology 
          // constraints
          LodNode* node = new LodNode(sequence, seqStart + tryPos, 
                                      seqStart + tryPos);
          nodeList->push_back(node);                    
        }
        colIt->toRight();
      } 
      while (colIt->lastColumn() == false && tryPos < maxTry);      
    }
  }
}

void LodGraph::createColumn(ColumnIteratorConstPtr colIt, NodeList* nodeList)
{

  // we don't add a new node for every palagous base in the reference genome
  // because we are not caching dupes within the iterator (not using
  // toRight() method).  So we only add a node for the current position

  const Sequence* refSequence = colIt->getReferenceSequence();
  
  hal_index_t pos = colIt->getReferenceSequencePosition();
  pos += refSequence->getStartPosition();

  LodNode* node = new LodNode(refSequence, pos, pos);
  _adjTable.addNode(node, colIt); 
  nodeList->push_back(node);  
}

void LodGraph::optimizeByExtension()
{
  for (GenomeNodesIterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    NodeList* nodeList = gni->second;

    for (NodeIterator i = nodeList->begin(); i != nodeList->end(); ++i)
    {
      LodNode* node = *i;
      node->extend(_extendFraction);
    }
  }
}

void LodGraph::optimizeByInsertion()
{
  vector<LodNode*> nodeBuffer;
  for (GenomeNodesIterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    NodeList* nodeList = gni->second;

    for (NodeIterator i = nodeList->begin(); i != nodeList->end(); ++i)
    {
      LodNode* node = *i;
      node->fillInEdges(nodeBuffer);      
    }
  }

  // these are all the nodes created by LodNode's fill-in method
  // we have to append them into the right litsts.
  // Note that these appends will result in the lists being unsorted.
  for (vector<LodNode*>::iterator i = nodeBuffer.begin(); 
       i != nodeBuffer.end(); ++i)
  {
    LodNode* newNode = *i;
    const Genome* genome = newNode->getSequence()->getGenome();
    GenomeNodesIterator gni = _genomeNodes.find(genome);
    assert(gni != _genomeNodes.end());
    NodeList* nodeList = gni->second;
    nodeList->push_back(newNode);
  }
}

// Note we assume nodelists are sorted
// Function should only be called in debug mode.
bool LodGraph::checkCoverage() const
{
  for (GenomeNodes::const_iterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    const Genome* genome = gni->first;    
    NodeList* nodeList = gni->second;
    hal_index_t prevPosition = -1;
    for (NodeIterator ni = nodeList->begin(); ni != nodeList->end(); ++ni)
    {
      LodNode* node = *ni;
      if (node->getSequence()->getGenome() != genome) 
      {
        cerr << "Node outside of genome " << genome->getName() << ": "
             << *node << endl;
        return false;
      }
      if (node->getStartPosition() != prevPosition + 1)
      {
        cerr << "Coverage gap: prev=" << prevPosition << ": "
             << *node << endl;
        return false;
      }      
      prevPosition = node->getEndPosition();
    }
    if ((hal_size_t)prevPosition != genome->getSequenceLength() - 1)
    {
      cerr << "Last node not at end.  Should be " <<
         genome->getSequenceLength() - 1 << ": "
           << prevPosition << endl;
      return false;
    }
  }
  return true;
}

void LodGraph::printDimensions(ostream& os) const
{
  for (GenomeNodes::const_iterator gni = _genomeNodes.begin(); 
       gni != _genomeNodes.end(); ++gni)
  {
    const Genome* genome = gni->first;    
    NodeList* nodeList = gni->second;
    os << genome->getName() << ": " << "nodeCount=" << nodeList->size() << " ";

    hal_size_t edgeCount = 0;
    hal_size_t maxDegree = 0;
    hal_size_t minDegree = numeric_limits<hal_size_t>::max();

    hal_size_t totalEdgeLength = 0;
    hal_size_t maxEdgeLength = 0;
    hal_size_t minEdgeLength = numeric_limits<hal_size_t>::max();

    hal_size_t fMin;
    hal_size_t fMax;
    hal_size_t fTot;
    hal_size_t rMin;
    hal_size_t rMax;
    hal_size_t rTot;
    
    for (NodeIterator ni = nodeList->begin(); ni != nodeList->end(); ++ni)
    {
      hal_size_t degree = (*ni)->getDegree();
      edgeCount += degree;
      maxDegree = max(degree, maxDegree);
      minDegree = min(degree, minDegree);
      assert(degree > 0);

      (*ni)->getEdgeLengthStats(fMin, fMax, fTot, rMin, rMax, rTot);
      totalEdgeLength += fTot + rTot;
      maxEdgeLength = max(maxEdgeLength, max(fMax, rMax));
      minEdgeLength = min(minEdgeLength, min(fMin, rMin));
    }    
    os << "edgeCount=" << edgeCount / 2 << " minDeg=" << minDegree << " "
       << "maxDeg=" << maxDegree << " totLen=" << totalEdgeLength << " "
       << "minLen=" << minEdgeLength << " maxLen=" << maxEdgeLength << endl;
  }
}
