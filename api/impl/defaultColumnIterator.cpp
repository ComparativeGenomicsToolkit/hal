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
:
  _root(root),
  _reference(reference),
  _index(columnIndex),
  _maxInsertionLength(maxInsertLength)
{
  init();
}
   
DefaultColumnIterator::~DefaultColumnIterator()
{
  
}

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

void DefaultColumnIterator::toRight() const
{
  resetColMap();
  assert(_activeStack.empty() == false);
  const Genome* refGenome = _activeStack.top();
  if (refGenome->getNumTopSegments() > 0)
  {
    TopMap::iterator i = _topMap.find(refGenome);
    assert(i != _topMap.end());
    LinkedTopIteratorPtr topIt = i->second;
    assert(topIt->_it.get() != NULL);
    hal_index_t position = topIt->_it->getStartPosition() + 1;
    topIt->_it->toRight(position);
    topIt->_dna->jumpTo(position);
    ColumnMap::iterator colIt = _colMap.find(refGenome);
    assert(colIt != _colMap.end());
    colIt->second.insert(topIt->_dna);
    // don't call recurisive functions if we've hit end of array
    // (todo: revise to handle insertions after reference)
    if ((hal_size_t)topIt->_it->getTopSegment()->getArrayIndex() <
        refGenome->getNumTopSegments())
    {
      updateParent(topIt);
      updateNextTopDup(topIt);
      updateParseDown(topIt);
    }
  } 
  else
  {
    assert(refGenome->getNumBottomSegments() > 0);
    BottomMap::iterator i = _bottomMap.find(refGenome);
    assert(i != _bottomMap.end());
    LinkedBottomIteratorPtr bottomIt = i->second;
    assert(bottomIt->_it.get() != NULL);
    hal_index_t position = bottomIt->_it->getStartPosition() + 1;
    bottomIt->_it->toRight(position);
    bottomIt->_dna->jumpTo(position);
    ColumnMap::iterator colIt = _colMap.find(refGenome);
    assert(colIt != _colMap.end());
    colIt->second.insert(bottomIt->_dna);
    hal_size_t numChildren = refGenome->getNumChildren();
    bottomIt->_children.resize(numChildren);
    // don't call recurisive functions if we've hit end of array
    // (todo: revise to handle insertions after reference)
    if ((hal_size_t)bottomIt->_it->getBottomSegment()->getArrayIndex() <
        refGenome->getNumBottomSegments())
    {
      for (size_t child = 0; child < numChildren; ++child)
      {
        updateChild(bottomIt, child);
      }
    }
  }  
  if (refGenome == _reference)
  {
    ++_index;
  }
}

// 
bool DefaultColumnIterator::leftOf(ColumnIteratorConstPtr other) const
{
  ColumnMap::const_iterator thisIt = _colMap.find(_reference);
  const ColumnMap* otherMap = other->getColumnMap();
  assert (_reference == other->getReferenceGenome());
  ColumnMap::const_iterator otherIt = otherMap->find(
    other->getReferenceGenome());
  if (thisIt != _colMap.end() && otherIt != otherMap->end() &&
      thisIt->second.empty() == false && otherIt->second.empty() == false)
  {
    DNASet::const_reverse_iterator thisDNA = thisIt->second.rbegin();
    DNASet::const_iterator otherDNA = otherIt->second.begin();
    return (*thisDNA)->leftOf(*otherDNA);
  }
  return false;
}


const Genome* DefaultColumnIterator::getReferenceGenome() const 
{
  return _reference;
}

const DefaultColumnIterator::ColumnMap* DefaultColumnIterator::getColumnMap() 
const
{
  return &_colMap;
}

void DefaultColumnIterator::init() const
{ 
  _colMap.clear();
  _topMap.clear();
  _bottomMap.clear();
  _topVisited.clear();
  _bottomVisited.clear();
  _activeStack = stack<const Genome*>();
  _activeStack.push(_reference);

  // if no root specified, we walk up to the alignment's root
  if (_root == NULL)
  {
    _root = _reference;
    while (_root->getParent() != NULL)
    {
      _root = _root->getParent();
    }
  }

  // Allocate a pair of linked iterators for every genome
  // todo: change!!! we don't need!! we just need iterators for 
  // genomes in the active stack!! 
  deque<const Genome*> bfQueue;
  bfQueue.push_back(_root);
  while (bfQueue.empty() == false)
  {
    const Genome* genome = bfQueue.front();
    bfQueue.pop_front();
    _topMap.insert(pair<const Genome*, LinkedTopIteratorPtr>(
                     genome, LinkedTopIteratorPtr(new LinkedTopIterator())));
    _bottomMap.insert(pair<const Genome*, LinkedBottomIteratorPtr>(
                        genome, 
                        LinkedBottomIteratorPtr(new LinkedBottomIterator())));
    _colMap.insert(pair<const Genome*, DNASet>(genome, DNASet()));
    hal_size_t numChildren = genome->getNumChildren();
    for (hal_size_t i = 0; i < numChildren; ++i)
    {
      bfQueue.push_back(genome->getChild(i));
    }
  }

  // dig out the reference genome and find the first iterators
  // move them to the start coordinate and recursively call the updates.  
  const Genome* refGenome = _activeStack.top();
  if (refGenome->getNumTopSegments() > 0)
  {
    TopMap::iterator i = _topMap.find(refGenome);
    assert(i != _topMap.end());
    LinkedTopIteratorPtr topIt = i->second;
    topIt->_it = refGenome->getTopSegmentIterator();
    topIt->_it->toSite(_index, true);
    topIt->_dna = refGenome->getDNAIterator(_index);
    ColumnMap::iterator colIt = _colMap.find(refGenome);
    assert(colIt != _colMap.end());
    colIt->second.insert(topIt->_dna);
    if (_index >= 0 && _index < (hal_index_t)refGenome->getSequenceLength())
    {
      updateParent(topIt);
      updateNextTopDup(topIt);
      updateParseDown(topIt);
    }
  } 
  else
  {
    assert(refGenome->getNumBottomSegments() > 0);
    BottomMap::iterator i = _bottomMap.find(refGenome);
    assert(i != _bottomMap.end());
    LinkedBottomIteratorPtr bottomIt = i->second;
    bottomIt->_it = refGenome->getBottomSegmentIterator();
    bottomIt->_it->toSite(_index, true);
    bottomIt->_dna = refGenome->getDNAIterator(_index);
    ColumnMap::iterator colIt = _colMap.find(refGenome);
    assert(colIt != _colMap.end());
    colIt->second.insert(bottomIt->_dna);
    hal_size_t numChildren = refGenome->getNumChildren();
    bottomIt->_children.resize(numChildren);
    if (_index >= 0 && _index < (hal_index_t)refGenome->getSequenceLength())
    {
      for (size_t child = 0; child < numChildren; ++child)
      {
        updateChild(bottomIt, child);
      }
    }
  }  
}

void DefaultColumnIterator::resetColMap() const
{
  for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i)
  {
    i->second.clear();
  }
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
    topIt->_parent->_dna->jumpTo( topIt->_parent->_it->getStartPosition());
    _colMap[parentGenome].insert(topIt->_parent->_dna);
    cout << "child parent " << topIt->_parent->_dna->getArrayIndex()
         << endl;

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
    _colMap[childGenome].insert(bottomIt->_children[index]->_dna);

    cout << "updating genome " << childGenome->getName() 
         << " (son of " << genome->getName() << ")"
         << " parent index " 
         << bottomIt->_dna->getArrayIndex()
         << " index " << bottomIt->_children[index]->_dna->getArrayIndex()
         << endl;
    
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
      
    _colMap[genome].insert(currentTopIt->_nextDup->_dna);

    // add the dup to the visit set if we're the reference
    // to do (add check for reference)
    const TopSegment* topSeg = currentTopIt->_nextDup->_it->getTopSegment();
    _topVisited.insert(VisitFlag(genome, topSeg->getArrayIndex()));
    
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
  // to do
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
    //_colMap[genome].push_back(bottomIt->_topParse->_dna);

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
    //_colMap[genome].push_back(topIt->_bottomParse->_dna);
    
    cout << "doing parse down on " << genome->getName()
         << " where incoming index is " << topIt->_dna->getArrayIndex()
         << " but outgoing index is " 
         << topIt->_bottomParse->_dna->getArrayIndex() << endl;

    // recurse on all the link's children
    for (hal_size_t i = 0; i <  topIt->_bottomParse->_children.size(); ++i)
    {
      updateChild(topIt->_bottomParse, i);
    }
  }
}  

