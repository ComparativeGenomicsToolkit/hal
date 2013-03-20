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
  for (EdgeIterator i = _edgeSet.begin(); i != _edgeSet.end(); ++i)
  {
    deleteEdge(i, false);
  }
  for (EdgeIterator i = _edgeSetZL.begin(); i != _edgeSetZL.end(); ++i)
  {
    deleteEdge(i, true);
  }
}

void LodNode::deleteEdge(EdgeIterator& edgeIt, bool zl)
{
  LodEdge* edge = *edgeIt;
  LodNode* other = edge->getOtherNode(this);
  EdgeIterator i = other->_edgeSet.find(edge);
  if (i != other->_edgeSet.end())
  {
    other->_edgeSet.erase(i);
    assert(other->_edgeSetZL.find(edge) == other->_edgeSetZL.end());
  }
  else
  {
    other->_edgeSetZL.erase(i);
    assert(other->_edgeSet.find(edge) == other->_edgeSet.end());
  }
  EdgeSet& edgeSet = zl ? _edgeSetZL : _edgeSet;
  edgeSet.erase(edgeIt);
  delete edge;
}

void addEdge(const Sequence* sequence, bool srcReversed, LodNode* tgt,
             bool tgtReversed)
{
  assert(_endPosition < tgt->_startPosition);
  hal_index_t length = tgt->_startPosition - _endPosition;
  assert(length >= 0);
  LodEdge* edge = new LodEdge(sequence, this, srcReversed, tgt,
                              tgtReversed);
  EdgeSet* srcSet;
  EdgeSet* tgtSet;
  if (length == 0)
  {
    srcSet = &_edgeSetZL;
    tgtSet = &tgt->_edgeSetZL;
  }
  else
  {
    srcSet = &_edgeSet;
    tgtSet = &tgt->_edgeSet;
  }
  assert(srcSet->find(edge) == srcSet->end() &&
         tgtSet->find(edge) == tgtSet->end());
  srcSet->insert(edge);
  tgtSet->insert(edge);
}
