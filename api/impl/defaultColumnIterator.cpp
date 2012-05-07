/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <deque>
#include "defaultColumnIterator.h"
#include "hal.h"

using namespace std;
using namespace hal;

DefaultColumnIterator::DefaultColumnIterator(const Genome* reference, 
                                             const Genome* root,
                                             hal_index_t columnIndex,
                                             hal_size_t maxInsertLength)
{

}
   
DefaultColumnIterator::~DefaultColumnIterator()
{

}

void DefaultColumnIterator::toRight() const
{

}

bool DefaultColumnIterator::equals(ColumnIteratorConstPtr other) const
{
  return false;
}


const Genome* DefaultColumnIterator::getReferenceGenome() const 
{
  return NULL;
}

void DefaultColumnIterator::init() const
{ 
  _topMap.clear();
  _bottomMap.clear();
  _topVisited.clear();
  _bottomVisited.clear();
  _activeStack = stack<const Genome*>();
  _activeStack.push(_reference);
  if (_root == NULL)
  {
    _root = _reference;
    while (_root->getParent() == NULL)
    {
      _root = _root->getParent();
    }
  }

  deque<const Genome*> bfQueue;
  bfQueue.push_back(_root);
  while (bfQueue.empty() == false)
  {
    const Genome* genome = bfQueue.front();
    bfQueue.pop_front();
    _topMap.insert(pair<const Genome*, LinkedTopIteratorPtr>(
                     genome, LinkedTopIteratorPtr()));
    _bottomMap.insert(pair<const Genome*, LinkedBottomIteratorPtr>(
                        genome, LinkedBottomIteratorPtr()));
    hal_size_t numChildren = genome->getNumChildren();
    for (hal_size_t i = 0; i < numChildren; ++i)
    {
      bfQueue.push_back(genome->getChild(i));
    }
  }
}

void DefaultColumnIterator::updateParent(LinkedTopIteratorPtr topIt) const
{
  if (topIt->_it->hasParent())
  {
    const Genome* genome = topIt->_it->getTopSegment()->getGenome();
    const Genome* parentGenome = genome->getParent();

    // no linked iterator for parent.  we create a new one and add the 
    // link in both directions
    if (topIt->_parent.get() == NULL)
    {
      assert(parentGenome != NULL);
      topIt->_parent = LinkedBottomIteratorPtr(new LinkedBottomIterator());
      topIt->_parent->_it = parentGenome->getBottomSegmentIterator();
      hal_size_t numChildren = parentGenome->getNumChildren();
      topIt->_parent->_children.resize(numChildren);
      for (hal_size_t i = 0; i < numChildren; ++i)
      {
        if (parentGenome->getChild(i) == genome)
        {
          topIt->_parent->_children[i] = topIt;
        } 
      }
    }

    // advance the parent's iterator to match topIt's (which should 
    // already have been updated. 
    topIt->_parent->_it->toParent(topIt->_it);

    // recurse on parent's parse edge
    updateParseUp(topIt->_parent);

    // recurse on parent's child edges (siblings to topIt)
    for (hal_size_t i = 0; i < topIt->_parent->_children.size(); ++i)
    {
      if (topIt->_parent->_children[i].get() == NULL ||
          topIt->_parent->_children[i]->_it->getTopSegment()->getGenome() != 
          genome)
      {
          updateChild(topIt->_parent, i);
      }
    }
  }  
}

void DefaultColumnIterator::updateChild(LinkedBottomIteratorPtr bottomIt, 
                                           hal_size_t index) const
{
  if (bottomIt->_it->hasChild(index))
  {
    assert(index < bottomIt->_children.size());
    const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();
    const Genome* childGenome = genome->getChild(index);

    // no linked iterator for child. we create a new one and add linke in
    // both directions 
    if (bottomIt->_children[index].get() == NULL)
    {
      assert(childGenome != NULL);
      bottomIt->_children[index] = LinkedTopIteratorPtr(new LinkedTopIterator());
      bottomIt->_children[index]->_it = childGenome->getTopSegmentIterator();
      bottomIt->_children[index]->_parent = bottomIt;
    }
    
    // advance the child's iterator to match bottomIt's (which should
    // have already been updated)
    bottomIt->_children[index]->_it->toChild(bottomIt->_it, index);

    //recurse on child's parse edge
    updateParseDown(bottomIt->_children[index]);
  }
}

void DefaultColumnIterator::updateNextTopDup(LinkedTopIteratorPtr topIt) const
{
  // to do
}

void DefaultColumnIterator::updateNextBottomDup(
  LinkedBottomIteratorPtr bottomIt) const
{
  // to do
}

void DefaultColumnIterator::updateParseUp(LinkedBottomIteratorPtr bottomIt)
   const
{
  if (bottomIt->_it->hasParseUp())
  {
    // no linked iterator for top parse, we create a new one
    if (bottomIt->_topParse.get() == NULL)
    {
      bottomIt->_topParse = LinkedTopIteratorPtr(new LinkedTopIterator());
      const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();
      bottomIt->_topParse->_it = genome->getTopSegmentIterator();
      bottomIt->_topParse->_bottomParse = bottomIt;
    }
    
    // advance the parse link's iterator to match bottomIt
    bottomIt->_topParse->_it->toParseUp(bottomIt->_it);

    // recurse on parse link's parent
    updateParent(bottomIt->_topParse);
  }
}
 
void DefaultColumnIterator::updateParseDown(LinkedTopIteratorPtr topIt) const
{
  if (topIt->_it->hasParseDown())
  {
    const Genome* genome = topIt->_it->getTopSegment()->getGenome();

    // no linked iterator for down parse, we create a new one
    if (topIt->_bottomParse.get() == NULL)
    {
      topIt->_bottomParse = LinkedBottomIteratorPtr(new LinkedBottomIterator());
      topIt->_bottomParse->_it = genome->getBottomSegmentIterator();
      topIt->_bottomParse->_topParse = topIt;
      hal_size_t numChildren = genome->getNumChildren();
      topIt->_bottomParse->_children.resize(numChildren);
      for (hal_size_t i = 0; i < numChildren; ++i)
      {
        if (genome->getChild(i) == genome)
        {
          topIt->_bottomParse->_children[i] = topIt;
        } 
      }
    }
    
    // advance the parse link's iterator to match topIt
    topIt->_bottomParse->_it->toParseDown(topIt->_it);
    
    // recurse on all the link's children
    for (hal_size_t i = 0; i <  topIt->_bottomParse->_children.size(); ++i)
    {
          updateChild(topIt->_bottomParse, i);
    }
  }
}  

