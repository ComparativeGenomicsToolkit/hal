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
#include "columnIteratorStack.h"
#include "hal.h"

using namespace std;
using namespace hal;

typedef ColumnIteratorStack::LinkedBottomIterator LinkedBottomIterator;
typedef ColumnIteratorStack::LinkedTopIterator LinkedTopIterator;
typedef ColumnIteratorStack::Entry Entry;

Entry::Entry(const Sequence* seq, hal_index_t first, hal_index_t index,
             hal_index_t last, hal_size_t size) : 
  _sequence(seq),
  _firstIndex(first),
  _index(index),
  _lastIndex(last),
  _cumulativeSize(size)
{
  _top._entry = this;
  _bottom._entry = this;
}

Entry::~Entry()
{
  freeLinks();
}

LinkedTopIterator* Entry::newTop()
{
  LinkedTopIterator* top = new LinkedTopIterator();
  top->_entry = this;
  _topLinks.push_back(top);
  return top;
}

LinkedBottomIterator* Entry::newBottom()
{
  LinkedBottomIterator* bottom = new LinkedBottomIterator();
  bottom->_entry = this;
  _bottomLinks.push_back(bottom);
  return bottom;
}

void Entry::freeLinks()
{
  size_t i;
  for (i = 0; i < _topLinks.size(); ++i)
  {
    delete _topLinks[i];
  }
  _topLinks.clear();
  _top._bottomParse = NULL;
  _top._parent = NULL;
  _top._nextDup = NULL;

  for (i = 0; i < _bottomLinks.size(); ++i)
  {
    delete _bottomLinks[i];
  }
  _bottomLinks.clear();
  _bottom._topParse = NULL;
  _bottom._children.clear();
}

ColumnIteratorStack::~ColumnIteratorStack()
{
  clear();
}

void ColumnIteratorStack::push(const Sequence* ref, hal_index_t index, 
                               hal_index_t lastIndex)
{
  assert(lastIndex >= index);
  assert(ref != NULL);
  hal_size_t cumulative = 0;
  if (_stack.size() > 0)
  {
    cumulative = top()->_cumulativeSize + lastIndex - index + 1;
  }
  Entry* entry = new Entry(ref, index, index, lastIndex, cumulative);
  _stack.push_back(entry);
}

void ColumnIteratorStack::pushStack(ColumnIteratorStack& otherStack)
{
  for (size_t i = 0; i < otherStack.size(); ++i)
  {
    _stack.push_back(otherStack[i]);
  }
  otherStack._stack.clear();
}

void ColumnIteratorStack::popDelete()
{
  assert(_stack.size() > 0);
  delete _stack.back();
  _stack.pop_back();
}

void ColumnIteratorStack::clear()
{
  while (_stack.size() > 0)
  {
    popDelete();
  }
}

void ColumnIteratorStack::resetLinks()
{
  for (size_t i = 0; i < _stack.size(); ++i)
  {
    _stack[i]->freeLinks();
  }
}
