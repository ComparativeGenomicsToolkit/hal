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


bool LodAdjTable::addNode(LodNode* node, ColumnIteratorConstPtr colIt)
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
  return true;
}

void LodAdjTable::writeAdjacenciesIntoNodes()
{
  for (Iterator i = _table.begin(); i != _table.end(); ++i)
  {
    const Sequence* sequence = i->first;
    RefSet* refSet = i->second;
    RefIterator cur = refSet->begin();
    RefIterator next = cur;
    RefIterator last = next;
    for (; cur != refSet->end(); ++cur)
    {
      // scan until we find a node with a bigger position
      for (; next->_pos == cur->_pos && next != refSet->end(); ++next);
   
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
        cur->_node->addEdge(sequence, cur->_reversed, last->_node, 
                            last->_reversed, distance - 1);
      }
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
