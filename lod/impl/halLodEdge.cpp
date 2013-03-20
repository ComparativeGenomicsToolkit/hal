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

LodEdge::LodEdge(const Sequnce* sequence, size_t length, LodNode* node1,
                 bool reversed1, LodNode* node2, bool reversed2) : 
  _sequence(sequence), _length(length), _node1(node1), _reversed1(reeversed1),
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
