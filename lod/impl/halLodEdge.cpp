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

LodEdge::LodEdge() : _sequence(NULL),
                     _node1(NULL), _pos1(NULL_INDEX), 
                     _node2(NULL), _pos2(NULL_INDEX), 
                     _reversed1(false), _reversed2(false)
{

}

LodEdge::LodEdge(const Sequence* sequence, 
                 LodNode* node1, hal_index_t pos1, bool reversed1,
                 LodNode* node2, hal_index_t pos2, bool reversed2) : 
  _sequence(sequence),
  _node1(node1), _pos1(pos1), 
  _node2(node2), _pos2(pos2), 
  _reversed1(reversed1), _reversed2(reversed2)
{
  assert(sequence != NULL);
  assert(_pos1 != _pos2);
  if (_pos1 > _pos2)
  {
    swap(_node1, _node2);
    swap(_pos1, _pos2);
    swap(_reversed1, _reversed2);
  }
  assert(!_node1 || _pos1 >= _sequence->getStartPosition());
  assert(!_node1 || _pos1 <= _sequence->getEndPosition());
  assert(!_node2 || _pos2 >= _sequence->getStartPosition());
  assert(!_node2 || _pos2 <= _sequence->getEndPosition());
  assert(_node1 || getLength() == 0);
  assert(_node2 || getLength() == 0);
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

void LodEdge::shrink(hal_size_t delta, bool left) 
{
  assert(delta <= getLength());
  if (left)
  {
    _pos1 += (hal_index_t)delta;
  }
  else
  {
    _pos2 -= (hal_index_t)delta;
  }
  assert(_pos1 < _pos2);
  assert(delta == 0 || (_node1 != NULL && _node2 != NULL));
  assert(!_node1 || _pos1 >= _sequence->getStartPosition());
  assert(!_node1 || _pos1 <= _sequence->getEndPosition());
  assert(!_node2 || _pos2 >= _sequence->getStartPosition());
  assert(!_node2 || _pos2 <= _sequence->getEndPosition());
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
  os << " endpoints=(" <<edge._pos1 << ", " << edge._pos2 
     << ") len=" << edge.getLength() 
     << " (" << edge._node1 << ", "
     << edge._reversed1 << ", " << edge._node2 << ", " << edge._reversed2
     << ")";
  return os;
}
