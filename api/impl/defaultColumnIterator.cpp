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

DefaultColumnIterator::DefaultColumnIterator(Genome* reference, 
                                             Genome* root,
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

void DefaultColumnIterator::toRightInParent(LinkedTopIteratorPtr topIt) const
{
  if (topIt->_parent->_it.get() != NULL)
  {
    LinkedBottomIteratorPtr& parent = topIt->_parent;
    // do torighton parent
  
    // parse up
    toRightInParseUp(parent);

    // descend down other children
    for (hal_size_t i = 0; i < parent->_children.size(); ++i)
    {
      if (parent->_children[i].get() != topIt.get())
      {
        toRightInChild(parent, i);
      }
    }  
  }
}

void DefaultColumnIterator::toRightInChild(LinkedBottomIteratorPtr bottomIt, 
                                           hal_size_t index) const
{
  assert(bottomIt->_children.size() > index);
  LinkedTopIteratorPtr& child = bottomIt->_children[index];
  // do toright on child

  // parse down
  toRightInParseDown(child);  
}

void DefaultColumnIterator::toRightInNextTopDup(LinkedTopIteratorPtr topIt) const
{
  // to do
}

void DefaultColumnIterator::toRightInNextBottomDup(
  LinkedBottomIteratorPtr bottomIt) const
{
  // to do
}

void DefaultColumnIterator::toRightInParseUp(LinkedBottomIteratorPtr bottomIt) const
{
  LinkedTopIteratorPtr& top = bottomIt->_topParse;
  // do toright in top

  toRightInParent(top);
}

void DefaultColumnIterator::toRightInParseDown(LinkedTopIteratorPtr topIt) const
{
  LinkedBottomIteratorPtr& bottom = topIt->_bottomParse;
  // do toright in bottom
  
  for (hal_size_t i = 0; i < bottom->_children.size(); ++i)
  {
    toRightInChild(bottom, i);
  }
}


