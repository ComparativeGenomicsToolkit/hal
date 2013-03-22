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

LodEdge::LodEdge() : _sequence(NULL), _start1(NULL), _length(0), _node1(NULL),
                     _node2(NULL), _reversed1(false), _reversed2(false),
                     _node1Left(NULL)
{

}

LodEdge::LodEdge(const Sequence* sequence, hal_index_t start1, bool node1Left,
                 size_t length, LodNode* node1,
                 bool reversed1, LodNode* node2, bool reversed2) : 
  _sequence(sequence), _start1(start1), _length(length), _node1(node1), 
  _node2(node2), _reversed1(reversed1), _reversed2(reversed2),
  _node1Left(node1Left)
{
  assert(sequence != NULL);
  assert(!_node1 || getStartPosition(_node1) >= _sequence->getStartPosition());
  assert(!_node2 || getStartPosition(_node2) >= _sequence->getStartPosition());
  assert(!_node1 || getStartPosition(_node1) <= _sequence->getEndPosition());
  assert(!_node2 || getStartPosition(_node2) <= _sequence->getEndPosition());
}

LodEdge::~LodEdge()
{
  
}

hal_index_t LodEdge::getStartPosition(const LodNode* node) const
{
  assert(node == _node1 || node == _node2);
  hal_index_t startPos = NULL_INDEX;
  if (node == _node1)
  {
    startPos = _start1;
  }
  else if (_node1Left == true)
  {
    startPos = _start1 + _length + 1;
  }
  else
  {
    startPos = _start1 - _length - 1;
  }
  assert(startPos >= _sequence->getStartPosition());
  assert(startPos <= _sequence->getEndPosition());
  return startPos;
}

LodNode* LodEdge::getOtherNode(const LodNode* node, bool* revThis,
                               bool* revOther)
{
  assert(node == _node1 || node == _node2);
  if (node == _node1)
  {
    if (revThis != NULL)
    {
      *revThis = _reversed1;
    }
    if (revOther != NULL)
    {
      *revOther = _reversed2;
    }
    return _node2;
  }
  else if (node == _node2)
  {
    if (revThis != NULL)
    {
      *revThis = _reversed2;
    }
    if (revOther != NULL)
    {
      *revOther = _reversed1;
    }
    return _node1;
  }
  return NULL;
}

void LodEdge::nullifyNode(const LodNode* node)
{
  assert(node == _node1 || node == _node2);
  if (node == _node1)
  {
    _node1 = NULL;
  }
  else
  {
    _node2 = NULL;
  }
}

ostream& hal::operator<<(ostream& os, const LodEdge& edge)
{
  os << "edge " << &edge << ": seq=";
  if (edge._sequence == NULL)
  {
    os << "NULL";
  }
  else
  {
    os << edge._sequence->getFullName();
  }
  os << " pos=" << edge._start1 << "(" << edge._node1Left 
     << ") len=" << edge.getLength() 
     << " (" << edge._node1 << ", "
     << edge._reversed1 << ", " << edge._node2 << ", " << edge._reversed2
     << ")";
  return os;
}

bool LodEdgePLess::operator()(const LodEdge* e1, const LodEdge* e2) const
{
  assert(e1 && e2);
  if (e1->_sequence < e2->_sequence)
  {
    return true;
  }
  else if (e1->_sequence == e2->_sequence)
  {
    if (e1->_node1 < e2->_node1)
    {
      return true;
    }
    else if (e1->_node1 == e2->_node1)
    {
      if (e1->_reversed1 == false && e2->_reversed1 == true)
      {
        return true;
      }
      else if (e1->_reversed1 == e2->_reversed1)
      {
        if (e1->_node2 < e2->_node2)
        {
          return true;
        }
        else if (e1->_node2 == e2->_node2)
        {
          return e1->_reversed2 == false && e2->_reversed2 == true;
        }
      }
    }
  }
  return false;
}
