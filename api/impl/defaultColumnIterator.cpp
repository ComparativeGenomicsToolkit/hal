/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <deque>
#include "defaultColumnIterator.h"
#include "hal.h"

using namespace std;
using namespace hal;

#ifdef _DISABLE_FOR_NOW
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

const ColumnIterator::SegmentMap& DefaultColumnIterator::getSegmentMap() const
{
  return _segmentMap;
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
  _activeStack.clear();
  _activeStack.push(_reference);
  if (_root == NULL)
  {
    _root = _referece;
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
    _topMap.insert(pair<const Genome*, LinkedTopIterator>(
                     genome, LinkedTopIterator()));
    _bottomMap.insert(pair<const Genome*, LinkedBottomIterator>(
                        genome, LinkedBottomIterator()));
    hal_size_t numChildren = genome->getNumChildren();
    for (hal_size_t i = 0; i < numChildren)
    {
      bfQueue.push_back(genome->getChild(i));
    }
  }
}

void DefaultColumnIterator::toRightInParent(LinkedTopIterator& topIt) const
{
  if (topIt._parent._it != NULL)
  {
    LinkedBottomIterator& parent - topIt._parent;
    // do torighton topIt._parent
  
    // parse up
    toRightInParseUp(parent);

    // descend down other children
    for (hal_size_t i = 0; i < parent._children.size(); ++i)
    {
      if (&parent._children[i] != & topIt)
      {
        toRightInChild(parent._children[i], i);
      }
    }  
  }
}

void DefaultColumnIterator::toRightInChild(LinkedBottomIterator& bottomIt, 
                                           hal_size_t index) const
{
  
}

void DefaultColumnIterator::toRightInNextTopDup(LinkedTopIterator& topIt) const
{

}

void DefaultColumnIterator::toRightInNextBottomDup(
  LinkedBottomIterator& bottomIt) const
{

}

void DefaultColumnIterator::toRightInParseUp(LinkedBottomIterator& bottomIt) const
{

}

void DefaultColumnIterator::toRightInParseDown(LinkedTopIterator& topIt) const
{

}

void DefaultColumnIterator::toRight() const
{
  deque<const Genome*> bfQueue;
  x
}
#endif
