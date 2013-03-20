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
    if (edge->getOtherNode(this) == NULL)
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
                      bool tgtReversed)
{
  assert(_endPosition < tgt->_startPosition);
  hal_index_t length = tgt->_startPosition - _endPosition;
  assert(length >= 0);
  LodEdge* edge = new LodEdge(sequence, length, this, srcReversed, tgt,
                              tgtReversed);
  _edges.push_back(edge);
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
  os << ")\n";

  hal_size_t ecount = 0;
  for (LodNode::EdgeList::const_iterator i = node._edges.begin(); 
       i != node._edges.end(); ++i)
  {
    os << ecount++ << ")" << *i << "; ";
  }
  os << "\n";
  return os;
}
