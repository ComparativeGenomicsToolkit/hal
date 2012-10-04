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
                                             hal_index_t lastColumnIndex,
                                             hal_size_t maxInsertLength,
                                             bool noDupes)
:
  _root(root),
  _maxInsertionLength(maxInsertLength),
  _noDupes(noDupes)
{
  assert (columnIndex >= 0 && lastColumnIndex >= columnIndex && 
          lastColumnIndex < (hal_index_t)reference->getSequenceLength());

  // allocate temp iterators
  _top = reference->getTopSegmentIterator(0);
  _next = _top->copy();

  // need to allocate the rearrangement from 
  if (reference->getParent() != NULL)
  {
    _rearrangement = reference->getRearrangement();
  }
  else if (reference->getNumChildren() > 0)
  {
    _rearrangement = reference->getChild(0)->getRearrangement();
  }
  _rearrangement->setAtomic(true);
  const Sequence* sequence = reference->getSequenceBySite(columnIndex);
  assert(sequence != NULL);

  // note columnIndex in genome (not sequence) coordinates
  _stack.push(sequence, columnIndex, lastColumnIndex);
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
  // compatible with old interface which allowed toRight() to go out
  // of bounds without crashing.
  if (_stack.size() == 1 && !_stack.topInBounds())
  {      
    return;
  }
  assert(_indelStack.size() == 0);
  
  do
  {
    // clean stack
    nextFreeIndex();
    while (_stack.size() > 1 && !_stack.topInBounds())
    {
      _stack.popDelete();
      nextFreeIndex();
    }

    // compatible with old interface which allowed toRight() to go out
    // of bounds without crashing.
    if (_stack.size() == 1 && !_stack.topInBounds())
    {      
      return;
    }

    _indelStack.clear();
    
    bool init = _stack.top()->_index == _stack.top()->_firstIndex ||
     (_stack.top()->_bottom._it.get() == NULL && 
      _stack.top()->_top._it.get() == NULL);

    recursiveUpdate(init);
    
    // move the index right
    ++_stack.top()->_index;

    // jump to next sequence in genome if necessary
    const Sequence* seq = _stack.top()->_sequence;
    if (_stack.size() == 1 && 
        _stack.top()->_index >= (hal_index_t)(seq->getStartPosition() + 
                                              seq->getSequenceLength()) &&
        _stack.top()->_index < (hal_index_t)(seq->getGenome()->getSequenceLength()))
    {
      _stack.top()->_sequence = 
         seq->getGenome()->getSequenceBySite(_stack.top()->_index);
      assert(_stack.top()->_sequence != NULL);
    }
  }
  while (_break == true);

  _ref = _stack.top()->_sequence;    

  // push the indel stack.  
  _stack.pushStack(_indelStack);

  // clean stack again
  nextFreeIndex();
  while (_stack.size() > 1 && !_stack.topInBounds())
  {
    _stack.popDelete();
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
     _stack.top()->_index > _stack.top()->_lastIndex;
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
  assert(_stack.size() > 0);
  return _stack.top()->_index;
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
  
  _stack.resetLinks();
}

// Starting from the reference sequence which is determined 
// from the stack, we start recursing over the entire column. 
// if init is specified, all the initial iterators are created
// then moved to the index (in the stack).  if init is false,
// all the existing iterators are moved to the right.
void DefaultColumnIterator::recursiveUpdate(bool init) const
{
/*  cout <<"update " << _stack.top()->_sequence->getName() << " "
       <<_stack.top()->_firstIndex << "," 
       <<_stack.top()->_index << ","
       <<_stack.top()->_lastIndex << endl;
*/

  resetColMap();
  _break = false;

  const Sequence* refSequence = _stack.top()->_sequence;
  const Genome* refGenome = refSequence->getGenome();
  if (refSequence->getNumTopSegments() > 0)
  {
    assert(_stack.size() > 0);
    LinkedTopIterator* topIt = &_stack.top()->_top;
    // first column, we search the genome for the site
    if (init == true)
    {    
      topIt->_it = refSequence->getTopSegmentIterator();
      topIt->_it->toSite(_stack.top()->_index, true);
      topIt->_dna = refGenome->getDNAIterator(_stack.top()->_index);
    }
    // otherwise, we scan forward from last visisted column
    else
    {
      assert(topIt->_it.get() != NULL);
      assert(topIt->_it->getReversed() == false);

      // catch up to nextfreeindex
      topIt->_it->slice(0, 0);
      while (topIt->_it->overlaps(_stack.top()->_index) == false)
      {
        topIt->_it->getReversed() ? topIt->_it->toLeft() : topIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - 
                                          topIt->_it->getStartPosition());
      topIt->_it->slice(offset, topIt->_it->getLength() - offset - 1);
      topIt->_dna->jumpTo(_stack.top()->_index);
    }
    assert(topIt->_it->getReversed() == false &&
           topIt->_dna->getReversed() == false);
    assert(topIt->_it->getStartPosition() == topIt->_dna->getArrayIndex());
    assert(topIt->_dna->getArrayIndex() == _stack.top()->_index);    
    assert(_stack.top()->_index <= _stack.top()->_lastIndex);
    assert(topIt->_it->getStartPosition() == topIt->_dna->getArrayIndex());

    if (colMapInsert(topIt->_dna) == false)
    {
      _break = true;
      return;
    }
    handleDeletion(topIt->_it);
    updateParent(topIt);
    updateNextTopDup(topIt);
    updateParseDown(topIt);
  } 

  else
  {
    assert(_stack.size() > 0);
    LinkedBottomIterator* bottomIt = &_stack.top()->_bottom;
    if (init == true)
    {
      bottomIt->_it = refSequence->getBottomSegmentIterator();
      bottomIt->_it->toSite(_stack.top()->_index, true);
      bottomIt->_dna = refGenome->getDNAIterator(_stack.top()->_index);
    }
    else
    {
      assert(bottomIt->_it.get() != NULL);
      assert(bottomIt->_it->getReversed() == false);

      // catch up to nextfreeindex
      bottomIt->_it->slice(0, 0);
      while (bottomIt->_it->overlaps(_stack.top()->_index) == false)
      {
        bottomIt->_it->getReversed() ? bottomIt->_it->toLeft() : 
           bottomIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - 
                                          bottomIt->_it->getStartPosition());
      bottomIt->_it->slice(offset, bottomIt->_it->getLength() - offset - 1);
      bottomIt->_dna->jumpTo(_stack.top()->_index);
    }

    assert(bottomIt->_it->getReversed() == false &&
           bottomIt->_dna->getReversed() == false);
    assert(bottomIt->_it->getStartPosition() == bottomIt->_dna->getArrayIndex());
    assert(bottomIt->_dna->getArrayIndex() == _stack.top()->_index);

    if (colMapInsert(bottomIt->_dna) == false)
    {
      _break = true;
      return;
    }
    hal_size_t numChildren = refSequence->getGenome()->getNumChildren();
    if (numChildren > bottomIt->_children.size())
    {
      bottomIt->_children.resize(numChildren, NULL);
    }
    assert(bottomIt->_it->getStartPosition() == 
           bottomIt->_dna->getArrayIndex());
    for (size_t child = 0; child < numChildren; ++child)
    {
      updateChild(bottomIt, child);
    }
  }
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
      const Genome* genome = _top->getTopSegment()->getGenome();
      const Genome* parent = genome->getParent();
      _top->slice(0, 0);
      assert(_rearrangement->getAtomic() == true);
      if (_rearrangement->identifyDeletionFromLeftBreakpoint(_top) == true && 
          _rearrangement->getLength() + _stack.top()->_cumulativeSize 
          <= _maxInsertionLength)
      {
        pair<hal_index_t, hal_index_t> deletedRange = 
           _rearrangement->getDeletedRange();
        assert((hal_size_t)(deletedRange.second - deletedRange.first) ==
               _rearrangement->getLength() - 1);

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
        _indelStack.push(bot->getBottomSegment()->getSequence(), 
                         deletedRange.first, deletedRange.second);

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
    _top->copy(inputTopIterator);
    bool reversed = _top->getReversed();
    // only handle an insertion if we are immediately left of the break       
    if (_top->getEndOffset() == 0)
    {
      _rearrangement->setAtomic(true);
      _top->slice(0, 0);
      _top->toRight();
      if (reversed == true)
      {
        _top->toReverse();
      }
      assert(_rearrangement->getAtomic() == true);
      if (_rearrangement->identifyInsertionFromLeftBreakpoint(_top) == true && 
          _rearrangement->getLength() + _stack.top()->_cumulativeSize 
          <= _maxInsertionLength)
      {
        pair<hal_index_t, hal_index_t> insertedRange = 
           _rearrangement->getInsertedRange();
        assert((hal_size_t)(insertedRange.second - insertedRange.first) ==
               _rearrangement->getLength() - 1);
        
/*
          cout << "\ninsertion found in " << inputTopIterator << endl;
          cout << "pushing " 
               << _top->getTopSegment()->getSequence()->getName()
               << "  " << insertedRange.first -
             _top->getTopSegment()->getSequence()->getStartPosition()<< " , " 
               << insertedRange.second -
             _top->getTopSegment()->getSequence()->getStartPosition()<< endl;
*/        
        _indelStack.push(_top->getTopSegment()->getSequence(), 
                         insertedRange.first, insertedRange.second);
      }
    }
  }
  return false;
}

void DefaultColumnIterator::updateParent(LinkedTopIterator* topIt) const
{
  const Genome* genome = topIt->_it->getTopSegment()->getGenome();
 
  if (!_break && genome != _root && topIt->_it->hasParent() && 
      checkRange(topIt->_dna))
  {
    const Genome* parentGenome = genome->getParent();

    // no linked iterator for parent.  we create a new one and add the 
    // link in both directions
    if (topIt->_parent == NULL)
    {
      assert(parentGenome != NULL);
      topIt->_parent = topIt->_entry->newBottom();
      topIt->_parent->_it = parentGenome->getBottomSegmentIterator();
      topIt->_parent->_dna = parentGenome->getDNAIterator();
      hal_size_t numChildren = parentGenome->getNumChildren();
      if (numChildren > topIt->_parent->_children.size())
      {
        topIt->_parent->_children.resize(numChildren, NULL);
      }
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
    if (colMapInsert(topIt->_parent->_dna) == false)
    {
      _break = true;
      return;
    }
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
      if (topIt->_parent->_children[i] == NULL ||
          topIt->_parent->_children[i]->_it->getTopSegment()->getGenome() != 
          genome)
      {
          updateChild(topIt->_parent, i);
      }
    }
  }  
}

void DefaultColumnIterator::updateChild(LinkedBottomIterator* bottomIt, 
                                        hal_size_t index) const
{
  if (!_break && bottomIt->_it->hasChild(index) && checkRange(bottomIt->_dna))
  {
    assert(index < bottomIt->_children.size());
    const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();
    const Genome* childGenome = genome->getChild(index);

    // no linked iterator for child. we create a new one and add linke in
    // both directions 
    if (bottomIt->_children[index] == NULL)
    {
      assert(childGenome != NULL);
      bottomIt->_children[index] = bottomIt->_entry->newTop();
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
    if(colMapInsert(bottomIt->_children[index]->_dna) == false)
    {
      _break = true;
      return;
    }
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

void DefaultColumnIterator::updateNextTopDup(LinkedTopIterator* topIt) const
{
  assert (topIt->_it.get() != NULL);
  if (_break || _noDupes == true ||
      topIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX ||
      !checkRange(topIt->_dna))
  {
    return;
  }

  hal_index_t firstIndex = topIt->_it->getTopSegment()->getArrayIndex();
  LinkedTopIterator* currentTopIt = topIt;
  const Genome* genome =  topIt->_it->getTopSegment()->getGenome();

  do  
  {
     // no linked iterator for paralog. we create a new one and add link
    if (currentTopIt->_nextDup == NULL)
    {
      currentTopIt->_nextDup = currentTopIt->_entry->newTop();
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
    if (colMapInsert(currentTopIt->_nextDup->_dna) == false)
    {
      _break = true;
      return;
    }
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

void DefaultColumnIterator::updateParseUp(LinkedBottomIterator* bottomIt)
   const
{
  if (!_break && bottomIt->_it->hasParseUp() && checkRange(bottomIt->_dna))
  {
    const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();

    // no linked iterator for top parse, we create a new one
    if (bottomIt->_topParse == NULL)
    {
      bottomIt->_topParse = bottomIt->_entry->newTop();
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
 
void DefaultColumnIterator::updateParseDown(LinkedTopIterator* topIt) const
{
  if (!_break && topIt->_it->hasParseDown() && checkRange(topIt->_dna))
  {
    const Genome* genome = topIt->_it->getTopSegment()->getGenome();

    // no linked iterator for down parse, we create a new one
    if (topIt->_bottomParse == NULL)
    {
      topIt->_bottomParse = topIt->_entry->newBottom();
      topIt->_bottomParse->_it = genome->getBottomSegmentIterator();
      topIt->_bottomParse->_dna = genome->getDNAIterator();
      topIt->_bottomParse->_topParse = topIt;
      hal_size_t numChildren = genome->getNumChildren();
      if (numChildren > topIt->_bottomParse->_children.size())
      {
        topIt->_bottomParse->_children.resize(numChildren, NULL);
      }
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

// moves index "right" until unvisited base is found
// if none exists in the current range, index is left one
// spot out of bounds (invalid) and return false.
void DefaultColumnIterator::nextFreeIndex() const
{
  hal_index_t index = _stack.top()->_index;

  VisitCache::iterator cacheIt = 
     _visitCache.find(_stack.top()->_sequence->getGenome());
  if (cacheIt != _visitCache.end())
  {
    PositionCache* posCache = cacheIt->second;
    bool found = posCache->find(index);
    while (found == true && index <= _stack.top()->_lastIndex)
    {
      ++index;
      found = posCache->find(index);
    }
  }
  _stack.top()->_index = index;
}

bool DefaultColumnIterator::colMapInsert(DNAIteratorConstPtr dnaIt) const
{
  const Sequence* sequence = dnaIt->getSequence();
  assert(sequence != NULL);
  
  // All reference bases need to get added to the cache
  bool updateCache = sequence == _stack[0]->_sequence;
  if (_maxInsertionLength == 0)
  {
    // Unless we don't do indels.  Here we just add reference elements
    // that are to right of the starting point
    assert (_stack.size() == 1);
    updateCache = _stack.top()->_index < dnaIt->getArrayIndex();
  }
  for (size_t i = 1; i < _stack.size() && !updateCache; ++i)
  {
    if (sequence == _stack[i]->_sequence)
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
      PositionCache* newSet = new PositionCache();
      cacheIt = _visitCache.insert(pair<const Genome*, PositionCache*>(
                                     sequence->getGenome(), newSet)).first;
    }
    found = cacheIt->second->insert(dnaIt->getArrayIndex()) == false;
  }
  else
  {
    found = cacheIt != _visitCache.end() && 
       cacheIt->second->find(dnaIt->getArrayIndex()) == true;
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

  return !found;
}

bool DefaultColumnIterator::checkRange(DNAIteratorConstPtr dnaIt) const
{
  StackEntry* entry = _stack.top();
  if (dnaIt->getSequence() == entry->_sequence)
  {
    assert (entry->_sequence->getGenome() == dnaIt->getGenome());
    // todo: check should be two-sided
    return dnaIt->getArrayIndex() >= entry->_index;
  }
  return true;
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
