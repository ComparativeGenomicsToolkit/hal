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
                                             hal_size_t maxInsertLength,
                                             bool endIterator)
:
  _root(root),
  _maxInsertionLength(maxInsertLength)
{
  // note columnIndex in genome (not sequence) coordinates
  init(reference, columnIndex, endIterator);
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
  recursiveUpdate(false);
}

// NOTE, will need changing once the stack gets used 
// ACTUALLY ITS BROKEN NOW 
bool DefaultColumnIterator::leftOf(ColumnIteratorConstPtr other) const
{
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

void DefaultColumnIterator::init(const Sequence* ref, hal_index_t index,
                                 bool endIterator) const
{ 
  eraseColMap();
  while (!_stack.empty())
  {
    _stack.pop();
  }

  // put ther reference on the stack
  StackEntry se = {ref,
                   index,
                   LinkedTopIteratorPtr(new LinkedTopIterator()),
                   LinkedBottomIteratorPtr(new LinkedBottomIterator()),
                   VisitSet()};
  _stack.push(se);

  if (endIterator == false)
  {  
    // might be an end interator.  in that case, do not recurse. 
    if ((hal_size_t)index < ref->getStartPosition() + ref->getSequenceLength())
    {
      recursiveUpdate(true);
    }
  }
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
  if (refSequence->getNumTopSegments() > 0)
  {
    assert(!_stack.empty());
    LinkedTopIteratorPtr topIt = _stack.top()._top;
    if (init == true)
    {    
      topIt->_it = refSequence->getTopSegmentIterator();
      topIt->_it->toSite(_stack.top()._index, true);
      topIt->_dna = refSequence->getDNAIterator(_stack.top()._index);
    }
    else
    {
      assert(topIt->_it.get() != NULL);
      assert(topIt->_it->getReversed() == false);
      _stack.top()._index = topIt->_it->getStartPosition() + 1;
      // do not handle iterating over multiple reference sequnces for now
      assert(topIt->_it->getTopSegment()->getSequence() == refSequence);
      topIt->_it->toRight(_stack.top()._index);
      topIt->_dna->jumpTo(_stack.top()._index);
    }
    colMapInsert(refSequence, topIt->_dna);
    if (_stack.top()._index >= 0 && 
        _stack.top()._index < (hal_index_t)refSequence->getSequenceLength())
    {
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
      bottomIt->_dna = refSequence->getDNAIterator(_stack.top()._index);
    }
    else
    {
      assert(bottomIt->_it.get() != NULL);
      assert(bottomIt->_it->getReversed() == false);
      _stack.top()._index = bottomIt->_it->getStartPosition() + 1;
      // do not handle iterating over multiple reference sequnces for now
      assert(bottomIt->_it->getBottomSegment()->getSequence() == refSequence);
      bottomIt->_it->toRight(_stack.top()._index);
      bottomIt->_dna->jumpTo(_stack.top()._index);
    }
    colMapInsert(refSequence, bottomIt->_dna);
    hal_size_t numChildren = refSequence->getGenome()->getNumChildren();
    bottomIt->_children.resize(numChildren);
    if (_stack.top()._index >= 0 && 
        _stack.top()._index < (hal_index_t)refSequence->getSequenceLength())
    {
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
 
  if (genome != _root && topIt->_it->hasParent())
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
/*
    cout << parentGenome->getName() << " " << topIt->_parent->_it->getReversed() << endl;

    cout << topIt->_parent->_it->getStartPosition() << " " << topIt->_it->getStartPosition() << endl << endl;
*/

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
  if (topIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX)
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

    // add the dup to the visit set if we're the reference
    // WORTHLESS, must do this every time! 
    const TopSegment* topSeg = currentTopIt->_nextDup->_it->getTopSegment();
    if (topSeg->getGenome() == _stack.top()._sequence->getGenome())
    {
      _stack.top()._visitSet.insert(topSeg->getArrayIndex());
    }
    
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
  if (bottomIt->_it->hasParseUp())
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

