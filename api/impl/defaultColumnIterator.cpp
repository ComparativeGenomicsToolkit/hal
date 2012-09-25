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
                                             bool endIterator)
:
  _root(root),
  _maxInsertionLength(maxInsertLength)
{
  // note columnIndex in genome (not sequence) coordinates
  init(reference, columnIndex, endIterator);

  // allocate temp iterators
  _top = reference->getTopSegmentIterator(0);
  _next = _top->copy();
}
   
DefaultColumnIterator::~DefaultColumnIterator()
{
  eraseColMap();
}

/*
hal_index_t DefaultColumnIterator::moveRightToNextUnvisited(
  LinkedTopIteratorPtr topIt) const
{
  // move one column to the right (as normal)
  hal_index_t position = topIt->_it->getStartPosition() + 1;
  topIt->_it->toRight(position);
  const TopSegment* topSeg = topIt->_it->getTopSegment();
  const Genome* genome = topSeg->getGenome();
  hal_size_t numSeg = genome->getNumTopSegments();

  VisitSet::iterator i;
  while((hal_size_t)topSeg->getArrayIndex() < numSeg)
  {  
    i = _topVisited.find(VisitFlag(genome, topSeg->getArrayIndex()));
    if (i != _topVisited.end())
    {
      // todo : we should be able to delete here (since we won't pass again)
      //_topVisited.remove(i);

      // move an entire segment to the right
      topIt->_it->slice(0, 0);
      topIt->_it->toRight();
    }
    else
    {
      break;
    }
  }
  
  if ((hal_size_t)topSeg->getArrayIndex() == numSeg)
  {
    return genome->getSequenceLength();
  }
  else
  {
    if (topIt->_it->getLength() > 1)
    {
      assert(topIt->_it->getStartOffset() == 0);
      topIt->_it->slice(0, topSeg->getLength() - 1);
    }
  }
  assert(topIt->_it->getLength() == 1);
  return topIt->_it->getStartPosition();
}
*/

void DefaultColumnIterator::toRight() const
{
  // we've hit the last index.  let's pop the stack. 
  // if we're at the bottom of the stack, something bad's happened
  // (maybe better to throw ane exception..)
  bool stackPushed = false;
  if (_stack.top()._index == _stack.top()._lastIndex && _stack.size() > 1)
  {
    _stack.pop();
    recursiveUpdate(true);
    cout << "popped\n\n";
  }
  else
  {
    cout << "calling HANDLE wwith " << stackPushed << " idx=" 
       << _stack.top()._index
       << endl;
    stackPushed = handleDeletion();
  }
  cout << "calling RU wwith " << stackPushed << " idx=" 
       << _stack.top()._index
       << endl;
  recursiveUpdate(stackPushed == true);
  cout << "after RU idx = " << _stack.top()._index << endl;

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
  return _stack.top()._sequence == _ref && 
     _stack.top()._index == _stack.top()._lastIndex;
}

// NOTE, will need changing once the stack gets used 
// ACTUALLY ITS BROKEN NOW 
bool DefaultColumnIterator::leftOf(ColumnIteratorConstPtr other) const
{
/*
  ColumnMap::const_iterator thisIt = _colMap.find(_stack.top()._sequence);
  const ColumnMap* otherMap = other->getColumnMap();
  assert (_stack.top()._sequence == other->getReferenceSequence());
  ColumnMap::const_iterator otherIt = otherMap->find(
    other->getReferenceSequence());
  if (thisIt != _colMap.end() && otherIt != otherMap->end() &&
      thisIt->second->empty() == false && otherIt->second->empty() == false)
  {
    DNASet::const_reverse_iterator thisDNA = thisIt->second->rbegin();
    DNASet::const_iterator otherDNA = otherIt->second->begin();
    return (*thisDNA)->leftOf(*otherDNA);
  }
*/
  return false;
}


const Genome* DefaultColumnIterator::getReferenceGenome() const 
{
  return _stack.top()._sequence->getGenome();
}

const Sequence* DefaultColumnIterator::getReferenceSequence() const 
{
  return _stack.top()._sequence;
}

const DefaultColumnIterator::ColumnMap* DefaultColumnIterator::getColumnMap() 
const
{
  return &_colMap;
}

hal_index_t DefaultColumnIterator::getArrayIndex() const
{
  assert(!_stack.empty());
  return _stack.top()._index;
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
}

void DefaultColumnIterator::pushStack(const Sequence* ref, 
                                      hal_index_t index,
                                      hal_index_t lastIndex,
                                      bool update) const
{
  // put ther reference on the stack
  RearrangementPtr rearrangement(NULL);

  // We only use rearrangements to detect deletions right now. Doing it in
  // both directions will really confuse things and doesn't fit the interface
  // as well.  This means that if the current reference sequence will only
  // have gaps in it related to ancestors.  If the reference is the root
  // then we can't even initialise the rearrangement object in the current
  // implementation.
  if (ref->getGenome()->getParent() != NULL)
  {
    rearrangement = ref->getRearrangement(0);
    rearrangement->setGapLengthThreshold(0);
  }

  StackEntry se = {ref,
                   index,
                   lastIndex,
                   false,
                   LinkedTopIteratorPtr(new LinkedTopIterator()),
                   LinkedBottomIteratorPtr(new LinkedBottomIterator()),
                   rearrangement,
                   VisitSet()};
  _stack.push(se);

  if (update == true)
  {  
    // might be an end interator.  in that case, do not recurse. 
    if ((hal_size_t)index < ref->getStartPosition() + ref->getSequenceLength())
    {
      recursiveUpdate(true);
    }
  }

  _ref= ref;
}

void DefaultColumnIterator::init(const Sequence* ref, hal_index_t index,
                                 bool endIterator) const
{ 
  eraseColMap();
  while (!_stack.empty())
  {
    _stack.pop();
  }

  const Genome* genome = ref->getGenome();
  hal_index_t lastIndex = (hal_index_t)genome->getSequenceLength() - 1;

  pushStack(ref, index, lastIndex, endIterator == false);
}

bool DefaultColumnIterator::handleDeletion() const
{
  if (_stack.top()._sequence->getGenome()->getParent() != NULL)
  {
    // test if we are the last base of a segment (
    // since there obviously can't be a deletion in any other
    // case:
    _top->copy(_stack.top()._top->_it);
    if (_top->getEndOffset() == 0)
    {
      const Genome* parent = _top->getTopSegment()->getGenome()->getParent();
      RearrangementPtr rearrangement = _stack.top()._rearrangement;
      if (_top->getReversed() == false)
      {
        cout << "testing " << _top << endl;
        _top->slice(0, 0);
        assert(rearrangement->getGapLengthThreshold() == 0);
        if (rearrangement->identifyDeletionFromLeftBreakpoint(_top) == true && 
            rearrangement->getLength() <= _maxInsertionLength)
        {
          BottomSegmentIteratorConstPtr bot = 
             parent->getBottomSegmentIterator(0);
          bot->toParent(_top);
          _next->copy(_top);
          if (_next->isLast() == false)
          {
            _next->toRight();
            hal_index_t start = _next->getStartPosition();
            hal_index_t last = start + _next->getLength() - 1;
            assert(last >+ start);
            cout << "pushing " << start << " , " << last << endl;
            pushStack(bot->getBottomSegment()->getSequence(), start, last, 
                      true);
            return true;
          }
        }
      }
      else
      {
        _top->copy(_stack.top()._top->_it);
        if (_top->isLast() == false)
        {
          // rearrangement doesn't work on inverted source sequence! 
          // so we do a test from the left neighbour
          _next->copy(_top);
          _next->toRight();
          _next->toReverse();
          if (rearrangement->identifyDeletionFromLeftBreakpoint(_next) == true 
              && rearrangement->getLength() <= _maxInsertionLength)
          {
            BottomSegmentIteratorConstPtr bot = 
               parent->getBottomSegmentIterator(0);
            bot->toParent(_top);
            hal_index_t start = bot->getStartPosition();
            hal_index_t last = _next->getStartPosition() - 1;
            assert(last < start);
            pushStack(bot->getBottomSegment()->getSequence(), start, last, 
                      true);
            return true;
          }
        }
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
//
// NOT CURRENTLY IMPLEMENTED:
//
// 1) Iterating across multiple sequences.
// 2) Properly detecting paralogies (!)
// 3) Insertions and stack updates
void DefaultColumnIterator::recursiveUpdate(bool init) const
{
  resetColMap();
  const Sequence* refSequence = _stack.top()._sequence;
  const Genome* refGenome = refSequence->getGenome();
  if (refSequence->getNumTopSegments() > 0)
  {
    assert(!_stack.empty());
    LinkedTopIteratorPtr topIt = _stack.top()._top;
    if (init == true)
    {    
      topIt->_it = refSequence->getTopSegmentIterator();
      topIt->_it->toSite(_stack.top()._index, true);
      topIt->_dna = refGenome->getDNAIterator(_stack.top()._index);
    }
    else
    {
      assert(topIt->_it.get() != NULL);
      assert(topIt->_it->getReversed() == false);
      nextFreeIndex();
      
      // do not handle iterating over multiple reference sequnces for now
      assert(topIt->_it->getTopSegment()->getSequence() == refSequence);
      if (_stack.top()._index < 0 || _stack.top()._index >= 
          (hal_index_t)(refSequence->getStartPosition() + 
                        refSequence->getSequenceLength()))
      {
        return;
      }
      topIt->_it->slice(0, 0);
      while (topIt->_it->overlaps(_stack.top()._index) == false)
      {
        topIt->_it->toRight();
      }
      hal_index_t offset = _stack.top()._index - topIt->_it->getStartPosition();
      assert(offset >= 0 && offset < (hal_index_t)topIt->_it->getLength());
      topIt->_it->slice(offset, topIt->_it->getLength() - offset - 1);
      topIt->_dna->jumpTo(_stack.top()._index);
    }
    colMapInsert(topIt->_dna, false);
    assert(topIt->_dna->getArrayIndex() == _stack.top()._index);
    
    if (_stack.top()._index >= 0 && 
        _stack.top()._index < (hal_index_t)(refSequence->getStartPosition() +
                                            refSequence->getSequenceLength()))
    {
      assert(topIt->_it->getStartPosition() == topIt->_dna->getArrayIndex());
      updateParent(topIt);
      updateNextTopDup(topIt);
      updateParseDown(topIt);
    }
  } 

  else
  {
    assert(!_stack.empty());
    LinkedBottomIteratorPtr bottomIt = _stack.top()._bottom;
    if (init == true)
    {
      bottomIt->_it = refSequence->getBottomSegmentIterator();
      bottomIt->_it->toSite(_stack.top()._index, true);
      bottomIt->_dna = refGenome->getDNAIterator(_stack.top()._index);
    }
    else
    {
      assert(bottomIt->_it.get() != NULL);
      assert(bottomIt->_it->getReversed() == false);
      nextFreeIndex();
      assert(bottomIt->_it->getBottomSegment()->getSequence() == refSequence);
      if (_stack.top()._index < 0 || _stack.top()._index >= 
          (hal_index_t)(refSequence->getStartPosition() + 
                        refSequence->getSequenceLength()))
      {
        return;
      }
      bottomIt->_it->slice(0, 0);
      while (bottomIt->_it->overlaps(_stack.top()._index) == false)
      {
        bottomIt->_it->toRight();
      }

      hal_index_t offset = _stack.top()._index -
         bottomIt->_it->getStartPosition();
      assert(offset >= 0 && offset < (hal_index_t)bottomIt->_it->getLength());
      bottomIt->_it->slice(offset, bottomIt->_it->getLength() - offset - 1);
      bottomIt->_dna->jumpTo(_stack.top()._index);
    }
    colMapInsert(bottomIt->_dna, false);
    assert(bottomIt->_dna->getArrayIndex() == _stack.top()._index);

    hal_size_t numChildren = refSequence->getGenome()->getNumChildren();
    bottomIt->_children.resize(numChildren);
    if (_stack.top()._index >= 0 && 
        _stack.top()._index <  (hal_index_t)(refSequence->getStartPosition() + 
                                             refSequence->getSequenceLength()))
    {
      assert(bottomIt->_it->getStartPosition() == 
             bottomIt->_dna->getArrayIndex());
      for (size_t child = 0; child < numChildren; ++child)
      {
        updateChild(bottomIt, child);
      }
    }
  }  
  
  // now that we moved from sets to vectors, we explictly sort.
  // (still on the fence on whether or not we need to keep sorted)
/*  for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i)
  {
    if (i->second->size() > 1)
    {
      sort(i->second->begin(), i->second->end());
    }
    }*/
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
/*
    cout << "updating genome " << childGenome->getName() 
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
  if (topIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX ||
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

void DefaultColumnIterator::updateNextBottomDup(
  LinkedBottomIteratorPtr bottomIt) const
{
  // don't think we need.  root never has any paralogies in current framweork. 
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

void DefaultColumnIterator::nextFreeIndex() const
{
  hal_index_t index = _stack.top()._index + 1;
  VisitSet::iterator i = _stack.top()._visitSet.find(index);
  while (i != _stack.top()._visitSet.end())
  {
//    _stack.top()._visitSet.erase(i);
    index += _stack.top()._reversed == false ? 1 : -1;
    i = _stack.top()._visitSet.find(index);
  }
  assert (index > _stack.top()._index);
  _stack.top()._index = index;
}

void DefaultColumnIterator::colMapInsert(DNAIteratorConstPtr dnaIt,
                                         bool updateVisitSet) const
{
  const Sequence* sequence = dnaIt->getSequence();
  assert(sequence != NULL);
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
  
  // add to the duplication visit set if it's the reference
  if (updateVisitSet == true &&
      sequence->getGenome() == _stack.top()._sequence->getGenome())
  {
    _stack.top()._visitSet.insert(dnaIt->getArrayIndex());
  }
}

bool DefaultColumnIterator::checkRange(DNAIteratorConstPtr dnaIt) const
{
  StackEntry& entry = _stack.top();
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
