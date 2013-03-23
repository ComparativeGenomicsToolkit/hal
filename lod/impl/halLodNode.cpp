/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cmath>
#include <cassert>
#include "halLodNode.h"
#include "halLodEdge.h"

using namespace std;
using namespace hal;

LodNode::LodNode() : _sequence(NULL), _startPosition(NULL_INDEX), 
                     _endPosition(NULL_INDEX)
{

}

LodNode::LodNode(const Sequence* sequence, hal_index_t start,
                 hal_index_t last) : _sequence(sequence), 
                                     _startPosition(start), 
                                     _endPosition(last)
{
  assert(_sequence != NULL);
  assert(_startPosition <= _endPosition);
  assert(_startPosition >= _sequence->getStartPosition());
  assert(_startPosition <= _sequence->getEndPosition());
  assert(_endPosition >= _sequence->getStartPosition());
  assert(_endPosition <= _sequence->getEndPosition());
}

LodNode::~LodNode()
{
  set<LodEdge*> tempSet;
  

  for (EdgeIterator i = _edges.begin(); i != _edges.end(); ++i)
  {    
    assert(tempSet.find(*i) == tempSet.end());
    tempSet.insert(*i);
    
    LodEdge* edge = *i;
    LodNode* other = edge->getOtherNode(this);
    if (other == NULL || other == this)
    {
      delete edge;
    }
    else
    {
      edge->nullifyNode(this);
    }
  }
}

void LodNode::addEdge(const Sequence* sequence,
                      hal_index_t srcPos, bool srcReversed, 
                      LodNode* tgt, hal_index_t tgtPos, bool tgtReversed)
{
  LodEdge* edge = new LodEdge(sequence, 
                              this, srcPos, srcReversed,
                              tgt, tgtPos, tgtReversed);
  assert(_edges.find(edge) == _edges.end());
  _edges.insert(edge);
  if (tgt != NULL && tgt != this)
  {
    assert(tgt->_edges.find(edge) == tgt->_edges.end());
    tgt->_edges.insert(edge);
  }
}

void LodNode::extend(double extendFraction)
{
  hal_size_t fMin;
  hal_size_t rMin; 
  getMinLengths(fMin, rMin);
  hal_size_t fLen = (hal_size_t)std::ceil((double)fMin * extendFraction);
  hal_size_t rLen = (hal_size_t)std::ceil((double)rMin * extendFraction);

  if (fLen > 0 || rLen > 0)
  {
    bool thisRev;
    LodNode* other;
    LodEdge* edge;
    for (EdgeIterator i = _edges.begin(); i != _edges.end(); ++i)
    {    
      edge = *i;
      other = edge->getOtherNode(this, &thisRev);
      bool left = edge->getStartPosition(this) < edge->getEndPosition(this);
      hal_size_t delta = thisRev ? rLen : fLen;

      edge->shrink(delta, left);
    }
    
    assert(_startPosition >= (hal_index_t)rLen);
    _startPosition -= (hal_index_t)rLen;
    
    assert(_endPosition + fLen <= (hal_size_t)_sequence->getEndPosition());
    _endPosition += (hal_index_t)fLen;
  }
}

void LodNode::fillInEdges(vector<LodNode*>& nodeBuffer)
{
  vector<LodEdge*> fBatch;
  vector<LodEdge*> rBatch;
  const Sequence* sequence = NULL;
  hal_index_t start = NULL_INDEX;
  EdgeIterator cur = _edges.begin();
  bool rev;
  while (true)
  {
    if (cur == _edges.end() || 
        (sequence != NULL && (*cur)->_sequence != sequence) ||
        (start != NULL_INDEX && (*cur)->getStartPosition(this) != start))
    {
      if (!fBatch.empty())
      {
        insertFillNodeIntoEdgeBatch(fBatch, nodeBuffer);
      }
      if (!rBatch.empty())
      {
        insertFillNodeIntoEdgeBatch(rBatch, nodeBuffer);
      }
      fBatch.clear();
      rBatch.clear();
    }
    if (cur == _edges.end())
    {
      break;
    }
    else if ((*cur)->getLength() > 0)
    {
      sequence = (*cur)->getSequence();
      start = (*cur)->getStartPosition(this);
      (*cur)->getOtherNode(this, &rev, NULL);
      if (rev)
      {
        assert(rBatch.empty() || 
               (*cur)->getLength() == rBatch.back()->getLength());
        rBatch.push_back(*cur);
      }
      else
      {
        assert(fBatch.empty() || 
               (*cur)->getLength()== fBatch.back()->getLength());
        fBatch.push_back(*cur);
      }
    }
    ++cur;
  }
}

void LodNode::insertFillNodeIntoEdgeBatch(vector<LodEdge*>& edgeBatch,
                                          vector<LodNode*>& nodeBuffer)
{
  assert(edgeBatch.size() > 0);
  LodEdge* edge = edgeBatch.at(0);
  assert(edge->getLength() > 0);
  hal_index_t firstStart = edge->getStartPosition(this);
  hal_index_t firstEnd = edge->getEndPosition(this);
  hal_index_t start = firstStart;
  hal_index_t end = firstEnd;
  if (start > end)
  {
    swap(start, end);
  }
  // change from exclusive edge coordinates to inclusive node coordinates
  start += 1;
  end -= 1;
  
  LodNode* newNode = new LodNode(edge->getSequence(), start, end);
  nodeBuffer.push_back(newNode);
  
  vector<LodEdge*>::iterator i;
  for (i = edgeBatch.begin(); i != edgeBatch.end(); ++i)
  {
    edge = *i;

    assert(edge->getStartPosition(this) == firstStart);
    assert(edge->getEndPosition((this)) == firstEnd);

    insertFillNodeIntoEdge(edge, newNode);
  }
}

void LodNode::insertFillNodeIntoEdge(LodEdge* edge, LodNode* newNode)
{
  // can do some of this in place and more efficiently.  start with
  // simplest way -- creating new edges from scratch -- to get it going.
  const Sequence* sequence = edge->getSequence();
  LodNode* leftNode = edge->_node1;
  LodNode* rightNode = edge->_node2;
  bool leftReversed = edge->_reversed1;
  bool rightReversed = edge->_reversed2;
  hal_index_t leftStart = edge->_pos1;
  hal_index_t rightStart  = edge->_pos2;
  
  // shouldn't be necessary but we do this so that future changesg 
  // to edge layout dont screw everything up
  if (leftStart > rightStart)
  {
    swap(leftNode, rightNode);
    swap(leftReversed, rightReversed);
    swap(leftStart, rightStart);
  }
  
  leftNode->_edges.erase(edge);
  rightNode->_edges.erase(edge);
  delete edge;
  
  leftNode->addEdge(sequence, leftStart, leftReversed,
                    newNode, leftStart + 1, !leftReversed);
  
  rightNode->addEdge(sequence, rightStart, rightReversed,
                     newNode, rightStart - 1, !rightReversed);
}

void LodNode::getEdgeLengthStats(hal_size_t& fMin, hal_size_t& fMax, 
                                 hal_size_t& fTot,
                                 hal_size_t& rMin, hal_size_t& rMax,
                                 hal_size_t& rTot) const
{
  fMin = numeric_limits<hal_size_t>::max();
  fMax = 0;
  fTot = 0;
  rMin = numeric_limits<hal_size_t>::max();
  rMax = 0;
  rTot = 0;  
  bool rev;

  for (EdgeConstIterator i = _edges.begin(); i != _edges.end(); ++i)
  {
    hal_size_t length = (*i)->getLength();
    (*i)->getOtherNode(this, &rev, NULL);
    if (rev == true)
    {
      rMin = min(length, rMin);
      rMax = max(length, rMax);
      rTot += length;
    }
    else
    {
      fMin = min(length, fMin);
      fMax = max(length, fMax);
      fTot += length;
    }
  }
}

void LodNode::getMinLengths(hal_size_t& fMin, hal_size_t& rMin) const
{
  fMin = numeric_limits<hal_size_t>::max();
  rMin = numeric_limits<hal_size_t>::max();  
  bool rev;

  for (EdgeConstIterator i = _edges.begin(); i != _edges.end(); ++i)
  {
    hal_size_t length = (*i)->getLength();
    (*i)->getOtherNode(this, &rev, NULL);
    if (rev == true)
    {
      rMin = min(length, rMin);
    }
    else
    {
      fMin = min(length, fMin);
    }
  }
}

ostream& hal::operator<<(ostream& os, const LodNode& node)
{
  os << "node " << &node << ": ";
  if (node.getSequence() != NULL)
  {
    os << node.getSequence()->getFullName();
  }
  else
  {
    os << "NULL";
  }
  os << "  (" << node.getStartPosition() << ", " << node.getEndPosition();
  os << ")";

  hal_size_t ecount = node.getDegree();
  if (ecount > 0)
  {
    os << "\n";
  }
  for (LodNode::EdgeList::const_iterator i = node._edges.begin(); 
       i != node._edges.end(); ++i)
  {
    LodNode* other = (*i)->getOtherNode(&node);
    if (other != NULL && other->getSequence() == node.getSequence())
    {
      os << " *";
    }
    else
    {
      os << "  ";
    }
    os << ecount++ << ")" << **i << "\n";
  }
  return os;
}
