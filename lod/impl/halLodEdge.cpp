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

LodEdge::LodEdge() : _sequence(NULL), _length(0), _node1(NULL),
                     _reversed1(false), _node2(NULL), _reversed2(false)
{

}

LodEdge::LodEdge(const Sequence* sequence, size_t length, LodNode* node1,
                 bool reversed1, LodNode* node2, bool reversed2) : 
  _sequence(sequence), _length(length), _node1(node1), _reversed1(reversed1),
  _node2(node2), _reversed2(reversed2)
{
  if (_node2 < _node1)
  {
    swap(_node1, _node2);
    swap(_reversed1, _reversed2);
  }
  assert(_node1 != _node2);
}

LodEdge::~LodEdge()
{
  
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
  os << "edge " << &edge << ":";
  if (edge.getSequence() != NULL)
  {
    os << edge.getSequence()->getFullName();
  }
  else
  {
    os << "NULL";
  }
  os << " len=" << edge.getLength() << " (" << edge._node1 << ", "
     << edge._reversed1 << ", " << edge._node2 << ", " << edge._reversed2
     << ")";
  return os;
}

/*LodEdgePLess::operator(const LodEdge* e1, const LodEdge* e2) const
{
  assert(e1 && e2);
  
  if (e1->_node1 < e2->_node1)
  {
    return true;
  }
  else if (e1->_node1 == e2->_node1)
  {
    if (e1->_node2 < e2->_node2)
    {
      return true;
    }
    else if (e1->_node2 == e2->_node2)
    {
      if (e1->_reversed1 == false && e2->_reversed1 == true)
      {
        return true;
      }
      else if (e1->_reversed1 == e2->_reversed1)
      {
        return e1->_reversed2 == false && e2->_reversed2 == true;
      }
    }
  }
  return false;
  }*/
