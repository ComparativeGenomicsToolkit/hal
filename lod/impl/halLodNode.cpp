/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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

void LodNode::addEdge(const Sequence* sequence, bool srcReversed, LodNode* tgt,
                      bool tgtReversed, hal_size_t length)
{
  assert(length >= 0);
  LodEdge* edge = new LodEdge(length, this, srcReversed, tgt,
                              tgtReversed);
  _edges.insert(edge);
  assert(tgt != this);
  tgt->_edges.insert(edge);
}

void LodNode::extend(bool reversed, hal_size_t length)
{
}

LodNode::EdgeIterator LodNode::getMinEdge(bool reversed)
{
  EdgeIterator minEdge = _edges.begin();
  EdgeIterator i = minEdge;
  bool rev;
  for (++i; i != _edges.end(); ++i)
  {
    (*i)->getOtherNode(this, &rev, NULL);
    if (rev == reversed && (*i)->getLength() < (*minEdge)->getLength())
    {
      minEdge = i;
    }
  }
  return minEdge;
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
