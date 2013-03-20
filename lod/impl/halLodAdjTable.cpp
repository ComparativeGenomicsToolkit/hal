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
      pair<Iterator, bool> res = _table.insert(sequence);
      if (res.second == true)
      {
        res.first->second = new refSset();
      }
      RefSet* refSet = res.first->second;
      ColumnIterator::DNASet::const_iterator dnaIt = dnaSet->begin();
      for (; dnaIt != dnaSet->end(); ++dnaIt)
      {
        NodeRef nodeRef(dnaIt->getArrayIndex(), dnaIt->getReversed(), node);
        insertRef(nodeRef, refSet);
      }      
    }
  }
}

void LodAdjTable::insertRef(const NodeRef& nodeRef, RefSet* refSet)
{
  pair<RefIterator, bool> res = refSet->insert(nodeRef);
  if (res.first == true)
  {
    bool overlaps = false;
    RefIterator i = res.first;
    if (i != refSet->begin())
    {
      --i;
      if (i->_pos + i->_node->getLength() >= nodeRef._pos)
      {
        overlaps = true;
      }
    }
    i = res.first;
    if (!overlaps && i != refSet->end())
    {
      ++i;
      if (i->_pos - i->_node->getLength() <= nodeRef._pos)
      {
        overlaps = true;
      }
    }
    if (overlaps == true)
    {
      refSet->erase(res.first);
    }
  }
}

void LodAdjTable::writeAdjacenciesIntoNodes()
{
  for (Iterator i = _table.begin(); i != _table.end(); ++i)
  {
    const Sequence* sequence = i->first;
    RefSet* refSet = i->second;
    for (RefIterator j = refSet->begin(); j != refSet->end(); ++j)
    {
      RefIterator next = j;
      ++next;
      if (next != refSet->end())
      {
        j->_node->addEdge(j->_reversed, next->_node, next->_reversed);
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
}
