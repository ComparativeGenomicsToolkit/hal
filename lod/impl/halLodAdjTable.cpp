/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include "halLodAdjTable.h"
#include "halLodNode.h"
#include "halLodEdge.h"
#include "halLodGraph.h"

using namespace std;
using namespace hal;

LodAdjTable::LodAdjTable()
{
  
}

LodAdjTable::~LodAdjTable()
{
  clear();
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
// TODO  
// 1) Test if total column size violates estimated
//    size bound.
// 2) Test if total similarity of added column to 
//    existing column is too great
// 3) Test if ref node hasn't already been added. 
////////////////////////////////////////////////////
////////////////////////////////////////////////////
bool LodAdjTable::canAddColumn(ColumnIteratorConstPtr colIt, hal_size_t step)
{
  return true;
}

void LodAdjTable::addNode(LodNode* node, ColumnIteratorConstPtr colIt)
{
  const ColumnIterator::ColumnMap* colMap = colIt->getColumnMap();
  ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
  for (; colMapIt != colMap->end(); ++colMapIt)
  {
    const ColumnIterator::DNASet* dnaSet = colMapIt->second;
    if (dnaSet && dnaSet->size() > 0)
    {
      const Sequence* sequence = colMapIt->first;
      pair<Iterator, bool> res = _table.insert(
        pair<const Sequence*, RefSet*>(sequence, NULL));
      if (res.second == true)
      {
        res.first->second = new RefSet();
      }
      RefSet* refSet = res.first->second;
      ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
      for (; dnaIt != dnaSet->end(); ++dnaIt)
      {
        NodeRef nodeRef((*dnaIt)->getArrayIndex(), (*dnaIt)->getReversed(), 
                        node);

        // Note: We assume nodes are length 1.  If they have been 
        // stretched, we probably need some sort of consistency check
        // to enforce that homologous nodes are the same size. 
        refSet->insert(nodeRef);
      }      
    }
  }
}

void LodAdjTable::writeAdjacenciesIntoNodes()
{
  for (Iterator i = _table.begin(); i != _table.end(); ++i)
  {
    RefSet* refSet = i->second;
    RefIterator cur = refSet->begin();
    RefIterator next;
    RefIterator last;
    for (; cur != refSet->end(); ++cur)
    {
      // scan until we find a node with a bigger position
      for (next = cur; next != refSet->end() && next->_pos == cur->_pos;
           ++next);
   
      // scan all nodes with this same position
      for (last = next; last != refSet->end() && last->_pos == next->_pos;
           ++last)
      {
        assert(last->_pos > cur->_pos);
        hal_size_t distance = last->_pos - cur->_pos;
        assert(distance > 0);
        
        // The length we pass is the number of spaces between the 
        // nodes.  The disatance we compute here is the difference
        // in their coordinates. 
        // last->_reversed is inverted because this node is to the right
        // and we are connecting to its "left side".  
        cur->_node->addEdge(cur->_reversed, last->_node, 
                            !last->_reversed, distance - 1);
      }
    }
    
    // put left caps (self loop on rev/rev with len 0)
    cur = refSet->begin();
    for (next = cur; next != refSet->end() && next->_pos == cur->_pos; ++next)
    {
      next->_node->addEdge(true, next->_node, true, 0);
    }
    
    // put right caps (self loop on for/for with len 0)
    RefRevIterator rcur = refSet->rbegin();
    for (RefRevIterator rnext = rcur; rnext != refSet->rend() && 
            rnext->_pos == rcur->_pos; ++rnext)
    {
      rnext->_node->addEdge(false, rnext->_node, false, 0);
    }
  }
}

void LodAdjTable::clear()
{
  for (Iterator i = _table.begin(); i != _table.end(); ++i)
  {
    delete i->second;
  }
  _table.clear();
}
