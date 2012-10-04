/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _COLUMNITERATORSTACK_H
#define _COLUMNITERATORSTACK_H

#include <set>
#include <stack>
#include <vector>
#include <map>
#include "halColumnIterator.h"
#include "halRearrangement.h"
#include "halCommon.h"

namespace hal {

class ColumnIteratorStack
{
public:

   struct Entry;
   struct LinkedTopIterator; 
   struct LinkedBottomIterator 
   {
      LinkedBottomIterator() : _topParse(NULL), _entry(NULL) {}
      BottomSegmentIteratorConstPtr _it;
      DNAIteratorConstPtr _dna;
      LinkedTopIterator* _topParse;
      std::vector<LinkedTopIterator*> _children;
      Entry* _entry;
   };

   struct LinkedTopIterator 
   {
      LinkedTopIterator() : _bottomParse(NULL), _parent(NULL), _nextDup(NULL),
                            _entry(NULL){}
      TopSegmentIteratorConstPtr _it;
      DNAIteratorConstPtr _dna;
      LinkedBottomIterator* _bottomParse;
      LinkedBottomIterator* _parent;
      LinkedTopIterator* _nextDup;
      Entry* _entry;
   };
   
   struct Entry 
   {
      Entry(const Sequence* seq, hal_index_t first, hal_index_t index,
                   hal_index_t last, hal_size_t size);
      ~Entry();
      LinkedTopIterator* newTop();
      LinkedBottomIterator* newBottom();
      void freeLinks();
      const Sequence* _sequence;
      hal_index_t _firstIndex;
      hal_index_t _index;
      hal_index_t _lastIndex;
      hal_size_t _cumulativeSize; 
      LinkedTopIterator _top;
      LinkedBottomIterator _bottom;
      std::vector<LinkedTopIterator*> _topLinks;
      std::vector<LinkedBottomIterator*> _bottomLinks;
   };
   
public:
   
   ~ColumnIteratorStack();
   void push(const Sequence* ref, hal_index_t index, hal_index_t lastIndex);
   void pushStack(ColumnIteratorStack& otherStack);
   void popDelete();
   void clear();
   Entry* top();
   Entry*& operator[](size_t i);
   size_t size() const;
   void resetLinks();
   bool topInBounds() const;

protected:

   std::vector<Entry*> _stack;
};

inline ColumnIteratorStack::Entry* ColumnIteratorStack::top()
{
  assert(_stack.size() > 0);
  return _stack.back();
}

inline ColumnIteratorStack::Entry*& ColumnIteratorStack::operator[](size_t i)
{
  assert(_stack.size() > i);
  return _stack[i];
}

inline size_t ColumnIteratorStack::size() const
{
  return _stack.size();
}

inline bool ColumnIteratorStack::topInBounds() const
{
  const Entry* e = _stack.back();
  return e->_index >= e->_firstIndex && e->_index <= e->_lastIndex;
}

}

#endif
