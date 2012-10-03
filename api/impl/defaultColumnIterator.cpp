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

DefaultColumnIterator::DefaultColumnIterator(const Sequence* reference, 
                                             const Genome* root,
                                             hal_index_t columnIndex,
                                             hal_index_t lastColumnIndex,
                                             hal_size_t maxInsertLength,
                                             bool noDupes)
:
  _root(root),
  _maxInsertionLength(maxInsertLength),
  _noDupes(noDupes)
{
  if (lastColumnIndex == NULL_INDEX)
  {
    lastColumnIndex = reference->getStartPosition() + 
       (hal_index_t)(reference->getSequenceLength() - 1);
  }

  // allocate temp iterators
  _top = reference->getTopSegmentIterator(0);
  _next = _top->copy();

  // note columnIndex in genome (not sequence) coordinates
  pushStack(_stack, reference, columnIndex, lastColumnIndex, false, true);
  toRight();
}
   
DefaultColumnIterator::~DefaultColumnIterator()
{
  eraseColMap();
  for (VisitCache::iterator i = _visitCache.begin();
       i != _visitCache.end(); ++i)
  {
    delete i->second;
  }
}

void DefaultColumnIterator::toRight() const
{
  // clean stack
  nextFreeIndex();
  while (_stack.size() > 1 && !inBounds())
  {
    _stack.pop_back();
    nextFreeIndex();
  }

  // compatible with old interface which allowed toRight() to go out
  // of bounds without crashing.
  if (_stack.size() == 1 && !inBounds())
  {      
    return;
  }

  assert(_indelStack.empty());
  
  size_t index = _stack.size() - 1;
  bool init = _stack.back()._index == _stack.back()._firstIndex ||
     (_stack.back()._bottom->_it.get() == NULL && 
      _stack.back()._top->_it.get() == NULL);

  recursiveUpdate(init);
  _ref = _stack.back()._sequence;
  
  // push the indel stack.  
  for (size_t i = 0; i < _indelStack.size(); ++i)
  {
    _stack.push_back(_indelStack[i]);
  }
  _indelStack.clear();

  // move the index right
  ++_stack[index]._index;

  // clean stack again
  nextFreeIndex();
  while (_stack.size() > 1 && !inBounds())
  {
    _stack.pop_back();
    nextFreeIndex();
  }

#ifndef NDEBUG
  set<pair<const Sequence*, hal_index_t> > coordSet;
  ColumnMap::const_iterator i, iNext;
  DNASet::const_iterator j;
  for (i = _colMap.begin(); i != _colMap.end(); ++i)
  {
    // check that the same coordinate not present for the same sequence
    for (j = i->second->begin(); j != i->second->end(); ++j)
    {
      pair<const Sequence*, hal_index_t> data(i->first, (*j)->getArrayIndex());
      assert(coordSet.insert(data).second == true);
    }
  }
#endif
}

bool DefaultColumnIterator::lastColumn() const
{
  return _stack.size() == 1 &&
     _stack.back()._index > _stack.back()._lastIndex;
}

const Genome* DefaultColumnIterator::getReferenceGenome() const 
{
  return _ref->getGenome();
}

const Sequence* DefaultColumnIterator::getReferenceSequence() const 
{
  return _ref;
}

const DefaultColumnIterator::ColumnMap* DefaultColumnIterator::getColumnMap() 
const
{
  return &_colMap;
}

hal_index_t DefaultColumnIterator::getArrayIndex() const
{
  assert(!_stack.empty());
  return _stack.back()._index;
}

void DefaultColumnIterator::defragment() const
{
  ColumnMap::iterator i = _colMap.begin();
  ColumnMap::iterator next;
  while ( i != _colMap.end())
  {
    next = i;
    ++next;
    if (i->second->empty())
    {
      delete i->second;
      _colMap.erase(i);
    }
    i = next;
  }
  
  ActiveStack::iterator j = _stack.begin();
  for (; j != _stack.end(); ++j)
  {
    j->_top->_bottomParse = LinkedBottomIteratorPtr();
    j->_top->_parent = LinkedBottomIteratorPtr();
    j->_top->_nextDup = LinkedTopIteratorPtr();
    j->_bottom->_topParse = LinkedTopIteratorPtr();
    j->_bottom->_children.clear();
    j->_bottom->_nextDup = LinkedBottomIteratorPtr();
  }
}

void DefaultColumnIterator::pushStack(ActiveStack& stack,
                                      const Sequence* ref, 
                                      hal_index_t index,
                                      hal_index_t lastIndex,
                                      bool reversed,
                                      bool update) const
{
  assert(lastIndex >= index);
  RearrangementPtr rearrangement(NULL);
  if (ref->getGenome()->getParent() != NULL)
  {
    rearrangement = ref->getRearrangement(0);
    rearrangement->setAtomic(true);
  }

  hal_size_t cumulativeSize = 0;
  if (_stack.size() > 0)
  {
    assert((hal_size_t)abs(lastIndex - index) <= _maxInsertionLength);
    cumulativeSize = _stack.back()._cumSize + abs(lastIndex - index + 1);
    assert(cumulativeSize <= _maxInsertionLength);
  }
  
//  cout << "pushing " << ref->getName() << " " <<index <<"," <<index <<","
//       <<lastIndex << endl;

  StackEntry se = {ref,
                   index,
                   index,
                   lastIndex,
                   cumulativeSize,
                   //reversed,
                   //DEUBUG
                   //
                   //
                   false,
                   LinkedTopIteratorPtr(new LinkedTopIterator()),
                   LinkedBottomIteratorPtr(new LinkedBottomIterator()),
                   rearrangement};
  stack.push_back(se);
}

bool DefaultColumnIterator::handleDeletion(TopSegmentIteratorConstPtr 
  inputTopIterator) const
{
  if (_maxInsertionLength > 0 && inputTopIterator->hasParent() == true)
  {
    _top->copy(inputTopIterator);
    bool reversed = _top->getReversed();
    if (reversed == true)
    {
      _top->toReverse();
    }
    // only handle a deletion if we are immediately left of the breakpoint   
    if (_top->getEndOffset() == 0)
    {
      const Genome* parent = _top->getTopSegment()->getGenome()->getParent();
      RearrangementPtr rearrangement = _stack.back()._rearrangement;
      _top->slice(0, 0);
      assert(rearrangement->getAtomic() == true);
      if (rearrangement->identifyDeletionFromLeftBreakpoint(_top) == true && 
          rearrangement->getLength() + _stack.back()._cumSize 
          <= _maxInsertionLength)
      {
        pair<hal_index_t, hal_index_t> deletedRange = 
           rearrangement->getDeletedRange();
        assert((hal_size_t)(deletedRange.second - deletedRange.first) ==
               rearrangement->getLength() - 1);

        BottomSegmentIteratorConstPtr bot = 
           parent->getBottomSegmentIterator(0);
        bot->toParent(_top);
        /*
        cout << "deletion found in " << bot << endl;
          cout << "pushing " 
               << bot->getBottomSegment()->getSequence()->getName()
               << "  " << deletedRange.first << " , " 
               << deletedRange.second << endl;
        */
          pushStack(_indelStack, bot->getBottomSegment()->getSequence(), 
                    deletedRange.first, deletedRange.second, 
                    reversed, true);

          return true;
      }
    }
  }
  return false;
}

bool DefaultColumnIterator::handleInsertion(TopSegmentIteratorConstPtr 
                                            inputTopIterator) const
{
  if (_maxInsertionLength > 0 && inputTopIterator->hasParent() == true)
  {
    // test if we are the last base of a segment (
    // since there obviously can't be a deletion in any other
    // case:
    _top->copy(inputTopIterator);
    bool reversed = _top->getReversed();
    if (reversed == true)
    {
      _top->toReverse();
    }
    // only handle an insertion if we are immediately left of the break       
    if (_top->getEndOffset() == 0)
    {
      RearrangementPtr rearrangement = _stack.back()._rearrangement;
      // test if next segment is inserted. if so, push to stack.
      _top->slice(0, 0);
      _top->toRight();
      assert(rearrangement->getAtomic() == true);
      if (rearrangement->identifyInsertionFromLeftBreakpoint(_top) == true && 
          rearrangement->getLength() + _stack.back()._cumSize 
          <= _maxInsertionLength)
      {
        pair<hal_index_t, hal_index_t> insertedRange = 
           rearrangement->getInsertedRange();
        assert((hal_size_t)(insertedRange.second - insertedRange.first) ==
               rearrangement->getLength() - 1);
        
/*
          cout << "\ninsertion found in " << inputTopIterator << endl;
          cout << "pushing " 
               << _top->getTopSegment()->getSequence()->getName()
               << "  " << insertedRange.first -
             _top->getTopSegment()->getSequence()->getStartPosition()<< " , " 
               << insertedRange.second -
             _top->getTopSegment()->getSequence()->getStartPosition()<< endl;
*/        
        pushStack(_indelStack, 
                  _top->getTopSegment()->getSequence(), insertedRange.first,
                  insertedRange.second, reversed, true);
      }
    }
  }
  return false;
}

// Starting from the reference sequence which is determined 
// from the stack, we start recursing over the entire column. 
// if init is specified, all the initial iterators are created
// then moved to the index (in the stack).  if init is false,
// all the existing iterators are moved to the right.
void DefaultColumnIterator::recursiveUpdate(bool init) const
{
/*  cout <<"update " << _stack.back()._sequence->getName() << " "
       <<_stack.back()._firstIndex << "," 
       <<_stack.back()._index << ","
       <<_stack.back()._lastIndex << endl;
*/

  resetColMap();

  const Sequence* refSequence = _stack.back()._sequence;
  const Genome* refGenome = refSequence->getGenome();
  if (refSequence->getNumTopSegments() > 0)
  {
    assert(!_stack.empty());
    LinkedTopIteratorPtr topIt = _stack.back()._top;
    // first column, we search the genome for the site
    if (init == true)
    {    
      topIt->_it = refSequence->getTopSegmentIterator();
      topIt->_it->toSite(_stack.back()._index, true);
      topIt->_dna = refGenome->getDNAIterator(_stack.back()._index);
      if (_stack.back()._reversed == true)
      {
        topIt->_it->toReverse();
        topIt->_dna->toReverse();
      }
    }
    // otherwise, we scan forward from last visisted column
    else
    {
      assert(topIt->_it.get() != NULL);
      assert(topIt->_it->getReversed() == _stack.back()._reversed);

      // catch up to nextfreeindex
      topIt->_it->slice(0, 0);
      while (topIt->_it->overlaps(_stack.back()._index) == false)
      {
        topIt->_it->getReversed() ? topIt->_it->toLeft() : topIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.back()._index - 
                                          topIt->_it->getStartPosition());
      topIt->_it->slice(offset, topIt->_it->getLength() - offset - 1);
      topIt->_dna->jumpTo(_stack.back()._index);
    }
    assert(topIt->_it->getReversed() == _stack.back()._reversed &&
           topIt->_dna->getReversed() == _stack.back()._reversed);
    assert(topIt->_it->getStartPosition() == topIt->_dna->getArrayIndex());
    assert(topIt->_dna->getArrayIndex() == _stack.back()._index);    
    assert(_stack.back()._index <= _stack.back()._lastIndex);
    assert(topIt->_it->getStartPosition() == topIt->_dna->getArrayIndex());

    colMapInsert(topIt->_dna);
    handleDeletion(topIt->_it);
    updateParent(topIt);
    updateNextTopDup(topIt);
    updateParseDown(topIt);
  } 

  else
  {
    assert(!_stack.empty());
    LinkedBottomIteratorPtr bottomIt = _stack.back()._bottom;
    if (init == true)
    {
      bottomIt->_it = refSequence->getBottomSegmentIterator();
      bottomIt->_it->toSite(_stack.back()._index, true);
      bottomIt->_dna = refGenome->getDNAIterator(_stack.back()._index);
      if (_stack.back()._reversed == true)
      {
        bottomIt->_it->toReverse();
        bottomIt->_dna->toReverse();
      }
    }
    else
    {
      assert(bottomIt->_it.get() != NULL);
      assert(bottomIt->_it->getReversed() == false);

      // catch up to nextfreeindex
      bottomIt->_it->slice(0, 0);
      while (bottomIt->_it->overlaps(_stack.back()._index) == false)
      {
        bottomIt->_it->getReversed() ? bottomIt->_it->toLeft() : 
           bottomIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.back()._index - 
                                          bottomIt->_it->getStartPosition());
      bottomIt->_it->slice(offset, bottomIt->_it->getLength() - offset - 1);
      bottomIt->_dna->jumpTo(_stack.back()._index);
    }

    assert(bottomIt->_it->getReversed() == _stack.back()._reversed &&
           bottomIt->_dna->getReversed() == _stack.back()._reversed);
    assert(bottomIt->_it->getStartPosition() == bottomIt->_dna->getArrayIndex());
    assert(bottomIt->_dna->getArrayIndex() == _stack.back()._index);

    colMapInsert(bottomIt->_dna);

    hal_size_t numChildren = refSequence->getGenome()->getNumChildren();
    bottomIt->_children.resize(numChildren);
    assert(bottomIt->_it->getStartPosition() == 
           bottomIt->_dna->getArrayIndex());
    for (size_t child = 0; child < numChildren; ++child)
    {
      updateChild(bottomIt, child);
    }
  }
}

void DefaultColumnIterator::resetColMap() const
{
  for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i)
  {
    i->second->clear();
  }
}

void DefaultColumnIterator::eraseColMap() const
{
  for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i)
  {
    delete i->second;
  }
  _colMap.clear();
}

void DefaultColumnIterator::updateParent(LinkedTopIteratorPtr topIt) const
{
  const Genome* genome = topIt->_it->getTopSegment()->getGenome();
 
  if (genome != _root && topIt->_it->hasParent() && 
      checkRange(topIt->_dna))
  {
    const Genome* parentGenome = genome->getParent();

    // no linked iterator for parent.  we create a new one and add the 
    // link in both directions
    if (topIt->_parent.get() == NULL)
    {
      assert(parentGenome != NULL);
      topIt->_parent = LinkedBottomIteratorPtr(new LinkedBottomIterator());
      topIt->_parent->_it = parentGenome->getBottomSegmentIterator();
      topIt->_parent->_dna = parentGenome->getDNAIterator();
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
    topIt->_parent->_dna->jumpTo( topIt->_parent->_it->getStartPosition());
    topIt->_parent->_dna->setReversed(topIt->_parent->_it->getReversed());
    colMapInsert(topIt->_parent->_dna);
    // cout << "child parent " << topIt->_parent->_dna->getArrayIndex() << endl;

    // recurse on parent's parse edge
    updateParseUp(topIt->_parent);
    if (topIt->_parent->_it->hasParseUp() &&
        topIt->_parent->_topParse->_it.get() != NULL)
    {
      handleDeletion(topIt->_parent->_topParse->_it);
    }

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
  if (bottomIt->_it->hasChild(index) && checkRange(bottomIt->_dna))
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
      bottomIt->_children[index]->_dna = childGenome->getDNAIterator();
      bottomIt->_children[index]->_parent = bottomIt;
    }
    
    // advance the child's iterator to match bottomIt's (which should
    // have already been updated)
    bottomIt->_children[index]->_it->toChild(bottomIt->_it, index);
    bottomIt->_children[index]->_dna->jumpTo(
      bottomIt->_children[index]->_it->getStartPosition());
    bottomIt->_children[index]->_dna->setReversed(
      bottomIt->_children[index]->_it->getReversed());
    colMapInsert(bottomIt->_children[index]->_dna);
    handleInsertion(bottomIt->_children[index]->_it);

/*    cout << "updating genome " << childGenome->getName() 
         << " (son of " << genome->getName() << ")"
         << " parent index " 
         << bottomIt->_dna->getArrayIndex()
         << " index " << bottomIt->_children[index]->_dna->getArrayIndex()
         << endl;
    cout << "parent it " << bottomIt->_it << endl;
    cout << "child it " << bottomIt->_children[index]->_it << endl << endl;
*/
    //recurse on paralgous siblings
    updateNextTopDup(bottomIt->_children[index]);

    //recurse on child's parse edge
    updateParseDown(bottomIt->_children[index]);
  }
}

void DefaultColumnIterator::updateNextTopDup(LinkedTopIteratorPtr topIt) const
{
  assert (topIt->_it.get() != NULL);
  if (_noDupes == true ||
      topIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX ||
      !checkRange(topIt->_dna))
  {
    return;
  }

  hal_index_t firstIndex = topIt->_it->getTopSegment()->getArrayIndex();
  LinkedTopIteratorPtr currentTopIt = topIt;
  const Genome* genome =  topIt->_it->getTopSegment()->getGenome();

  do  
  {
     // no linked iterator for paralog. we create a new one and add link
    if (currentTopIt->_nextDup.get() == NULL)
    {
      currentTopIt->_nextDup = LinkedTopIteratorPtr(new LinkedTopIterator());
      currentTopIt->_nextDup->_it = genome->getTopSegmentIterator();
      currentTopIt->_nextDup->_dna = genome->getDNAIterator();
      currentTopIt->_nextDup->_parent = currentTopIt->_parent;
    }
    
    // advance the dups's iterator to match currentTopIt's (which should
    // have already been updated)
    currentTopIt->_nextDup->_it = currentTopIt->_it->copy();
    currentTopIt->_nextDup->_it->toNextParalogy();
    currentTopIt->_nextDup->_dna->jumpTo(
      currentTopIt->_nextDup->_it->getStartPosition());
    currentTopIt->_nextDup->_dna->setReversed(
      currentTopIt->_nextDup->_it->getReversed());
    colMapInsert(currentTopIt->_nextDup->_dna);
    handleInsertion(currentTopIt->_nextDup->_it);
    
    // recurse on duplicate's parse edge
    updateParseDown(currentTopIt->_nextDup);

    // advance current it to the next paralog
    currentTopIt = currentTopIt->_nextDup;
  } 
  while (currentTopIt->_it->getTopSegment()->getNextParalogyIndex() != 
         NULL_INDEX &&
         currentTopIt->_it->getTopSegment()->getNextParalogyIndex() != 
         firstIndex);
}

void DefaultColumnIterator::updateParseUp(LinkedBottomIteratorPtr bottomIt)
   const
{
  if (bottomIt->_it->hasParseUp() && checkRange(bottomIt->_dna))
  {
    const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();

    // no linked iterator for top parse, we create a new one
    if (bottomIt->_topParse.get() == NULL)
    {
      bottomIt->_topParse = LinkedTopIteratorPtr(new LinkedTopIterator());
      bottomIt->_topParse->_it = genome->getTopSegmentIterator();
      bottomIt->_topParse->_dna = genome->getDNAIterator();
      bottomIt->_topParse->_bottomParse = bottomIt;
    }
    
    // advance the parse link's iterator to match bottomIt
    bottomIt->_topParse->_it->toParseUp(bottomIt->_it);
    bottomIt->_topParse->_dna->jumpTo(
      bottomIt->_topParse->_it->getStartPosition());
    bottomIt->_topParse->_dna->setReversed(
      bottomIt->_topParse->_it->getReversed());
    assert(bottomIt->_topParse->_dna->getArrayIndex() ==
           bottomIt->_dna->getArrayIndex());

    // recurse on parse link's parent
    updateParent(bottomIt->_topParse);

    //recurse on parse link's paralogous siblings
    updateNextTopDup(bottomIt->_topParse);
  }
}
 
void DefaultColumnIterator::updateParseDown(LinkedTopIteratorPtr topIt) const
{
  if (topIt->_it->hasParseDown() && checkRange(topIt->_dna))
  {
    const Genome* genome = topIt->_it->getTopSegment()->getGenome();

    // no linked iterator for down parse, we create a new one
    if (topIt->_bottomParse.get() == NULL)
    {
      topIt->_bottomParse = LinkedBottomIteratorPtr(new LinkedBottomIterator());
      topIt->_bottomParse->_it = genome->getBottomSegmentIterator();
      topIt->_bottomParse->_dna = genome->getDNAIterator();
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

    
    topIt->_bottomParse->_dna->jumpTo(
      topIt->_bottomParse->_it->getStartPosition());
    topIt->_bottomParse->_dna->setReversed(
      topIt->_bottomParse->_it->getReversed());
    assert(topIt->_bottomParse->_dna->getArrayIndex() ==
           topIt->_dna->getArrayIndex());

/*
    cout << "doing parse down on " << genome->getName()
         << " where incoming index is " << topIt->_dna->getArrayIndex()
         << " but outgoing index is " 
         << topIt->_bottomParse->_dna->getArrayIndex() << endl;
*/
    // recurse on all the link's children
    for (hal_size_t i = 0; i <  topIt->_bottomParse->_children.size(); ++i)
    {
      updateChild(topIt->_bottomParse, i);
    }
  }
}  


bool DefaultColumnIterator::inBounds() const
{
  StackEntry& e = _stack.back();
  return e._index >= e._firstIndex && e._index <= e._lastIndex;
}

// moves index "right" until unvisited base is found
// if none exists in the current range, index is left one
// spot out of bounds (invalid) and return false.
bool DefaultColumnIterator::nextFreeIndex() const
{
  hal_index_t index = _stack.back()._index;
  bool success = true;

  VisitCache::iterator cacheIt = 
     _visitCache.find(_stack.back()._sequence->getGenome());
  if (cacheIt != _visitCache.end())
  {
    VisitSet* visitSet = cacheIt->second;
    VisitSet::iterator i = visitSet->find(index);
    while (i != visitSet->end())
    {
      ++index;
      i = visitSet->find(index);
      if (index - _stack.back()._lastIndex == 1)
      {
        success = false;
        break;
      }
    }
  }
  _stack.back()._index = index;
  return success;
}

bool DefaultColumnIterator::colMapInsert(DNAIteratorConstPtr dnaIt) const
{
  const Sequence* sequence = dnaIt->getSequence();
  assert(sequence != NULL);
  
  bool updateCache = false;
  // Update cache Case 1: on base reference ahead not index
  // ie, we looped back onto the reference somehow via paralogy
  if (sequence->getGenome() == _stack[0]._sequence->getGenome() &&
       dnaIt->getArrayIndex() != _stack[0]._index)
  {
    updateCache = true;
  }
  // Update cache Case 2: the current base is in an indel.  (would be
  // nice to only check top of stack here but we do the whole thing for now)
  for (size_t i = 1; i < _stack.size() && !updateCache; ++i)
  {
    if (sequence == _stack[i]._sequence)
    {
      updateCache = true;
    }
  }

  bool found = false;
  VisitCache::iterator cacheIt = _visitCache.find(sequence->getGenome());
  if (updateCache == true)
  {
    if (cacheIt == _visitCache.end())
    {
      VisitSet* newSet = new VisitSet();
      cacheIt = _visitCache.insert(
        pair<const Genome*, VisitSet*>(sequence->getGenome(), newSet)).first;
    }
    found = cacheIt->second->insert(dnaIt->getArrayIndex()).second == false;
  }
  else
  {
    found = cacheIt != _visitCache.end() && 
       cacheIt->second->find(dnaIt->getArrayIndex()) != cacheIt->second->end();
  }

  // insert into the column data structure to pass out to client
  if (found == false)
  {
    ColumnMap::iterator i = _colMap.lower_bound(sequence);
    if(i != _colMap.end() && !(_colMap.key_comp()(sequence, i->first)))
    {
      i->second->push_back(dnaIt);
    }
    else
    {
      DNASet* dnaSet = new DNASet();
      dnaSet->push_back(dnaIt);
      _colMap.insert(i, ColumnMap::value_type(sequence, dnaSet));
    }
  }
  assert(found == false);

  return !found;
}

bool DefaultColumnIterator::checkRange(DNAIteratorConstPtr dnaIt) const
{
  StackEntry& entry = _stack.back();
  if (dnaIt->getSequence() == entry._sequence)
  {
    assert (entry._sequence->getGenome() == dnaIt->getGenome());
    // todo: check should be two-sided
    if (entry._reversed == false)
    {
      return dnaIt->getArrayIndex() >= entry._index;
    }
    else
    {
      return dnaIt->getArrayIndex() <= entry._index;
    }
  }
  return true;
}

