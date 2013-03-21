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
  assert((_startPosition == NULL_INDEX && _endPosition == NULL_INDEX) ||
         (_startPosition <= _endPosition && _sequence != NULL));
}

LodNode::~LodNode()
{
  for (EdgeIterator i = _edges.begin(); i != _edges.end(); ++i)
  {    
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

void LodNode::addEdge(bool srcReversed, LodNode* tgt,
                      bool tgtReversed, hal_size_t length)
{
  assert(length >= 0);
  LodEdge* edge = new LodEdge(length, this, srcReversed, tgt,
                              tgtReversed);
  _edges.insert(edge);
  tgt->_edges.insert(edge);
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
      if (thisRev == true)
      {
        edge->_length -= rLen;
      }
      else
      {
        edge->_length -= fLen;
      }
    }
    assert(_startPosition >= (hal_index_t)rLen);
    _startPosition -= (hal_index_t)rLen;
    
    assert(_endPosition + fLen < 
           _sequence->getStartPosition() + _sequence->getSequenceLength());
    _endPosition += (hal_index_t)fLen;
  }
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

  for (EdgeIterator i = _edges.begin(); i != _edges.end(); ++i)
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

  for (EdgeIterator i = _edges.begin(); i != _edges.end(); ++i)
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
    os << "  " << ecount++ << ")" << **i << "\n";
  }
  return os;
}
