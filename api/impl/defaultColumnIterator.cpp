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
                                             const set<const Genome*>* targets,
                                             hal_index_t columnIndex,
                                             hal_index_t lastColumnIndex,
                                             hal_size_t maxInsertLength,
                                             bool noDupes,
                                             bool noAncestors,
                                             bool reverseStrand,
                                             bool unique,
                                             bool onlyOrthologs)
:
  _maxInsertionLength(maxInsertLength),
  _noDupes(noDupes),
  _noAncestors(noAncestors),
  _reversed(reverseStrand),
  _tree(NULL),
  _unique(unique),
  _onlyOrthologs(onlyOrthologs)
{
  assert (columnIndex >= 0 && lastColumnIndex >= columnIndex && 
          lastColumnIndex < (hal_index_t)reference->getSequenceLength());

  // allocate temp iterators
  if (reference->getNumTopSegments() > 0)
  {
    _top = reference->getTopSegmentIterator(0);
    _next = _top->copy();
  }
  else if (reference->getChild(0) != NULL)
  {
    _top = reference->getChild(0)->getTopSegmentIterator(0);
    _next = _top->copy();
  }

  if (_maxInsertionLength > 0)
  {
    // need to allocate the rearrangement from 
    if (reference->getParent() != NULL)
    {
      _rearrangement = reference->getRearrangement(0, 0, 1., true);
    }
    else if (reference->getNumChildren() > 0)
    {
      _rearrangement = reference->getChild(0)->getRearrangement(0, 0, 1., true);
    }
  }
  const Sequence* sequence = reference->getSequenceBySite(columnIndex);
  assert(sequence != NULL);
  _ref =sequence;    

  // compute all genomes in search scope (spanning tree of reference and targets)
  // if targets is empty we just visit everything. 
  if (targets != NULL && !targets->empty())
  {
    _targets = *targets;
    _targets.insert(reference);
    getGenomesInSpanningTree(_targets, _scope);
  }
  
  // note columnIndex in genome (not sequence) coordinates
  _stack.push(sequence, columnIndex, lastColumnIndex);

  toRight();
}
   
DefaultColumnIterator::~DefaultColumnIterator()
{
  eraseColMap();
  clearVisitCache();
  clearTree();
}

void DefaultColumnIterator::toRight() const
{
  clearTree();

  // keep the current position so that when client calls
  // getReferenceXXX() methods, they get the state before 
  // toRight is called. 
  _prevRefSequence = _ref;
  _prevRefIndex = _stack[0]->_index - _ref->getStartPosition();

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
      _ref = _stack.top()->_sequence;    
    }
  }
  while (_break == true);

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

void DefaultColumnIterator::toSite(hal_index_t columnIndex, 
                                   hal_index_t lastColumnIndex,
                                   bool clearCache) const
{
  clearTree();

  const Genome* reference = getReferenceGenome();
  assert (columnIndex >= 0 && lastColumnIndex >= columnIndex && 
          lastColumnIndex < (hal_index_t)reference->getSequenceLength());  

  const Sequence* sequence = reference->getSequenceBySite(columnIndex);
  assert(sequence != NULL);
  _ref =sequence;    
  _stack.clear();
  _indelStack.clear();
  if (clearCache == true)
  {
      clearVisitCache();
  }
  defragment();
  // note columnIndex in genome (not sequence) coordinates
  _stack.push(sequence, columnIndex, lastColumnIndex);
  toRight();
  assert(getReferenceSequencePosition() + sequence->getStartPosition() == 
         columnIndex);
}

bool DefaultColumnIterator::lastColumn() const
{
  return _stack.size() == 1 &&
     _stack.top()->_index > _stack.top()->_lastIndex;
}

const Genome* DefaultColumnIterator::getReferenceGenome() const 
{
  return _prevRefSequence->getGenome();
}

const Sequence* DefaultColumnIterator::getReferenceSequence() const 
{
  return _prevRefSequence;
}

hal_index_t DefaultColumnIterator::getReferenceSequencePosition() const 
{
  return _prevRefIndex;
}

const DefaultColumnIterator::ColumnMap* DefaultColumnIterator::getColumnMap() 
const
{
  return &_colMap;
}

hal_index_t DefaultColumnIterator::getArrayIndex() const
{
  assert(_stack.size() > 0);
  return _stack[0]->_index;
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

bool DefaultColumnIterator::isCanonicalOnRef() const
{
  assert(_stack.size() > 0);
  assert(_leftmostRefPos >= 0 && (hal_size_t)_leftmostRefPos < 
         _stack[0]->_sequence->getGenome()->getSequenceLength());
  return _leftmostRefPos >= _stack[0]->_firstIndex &&
     _leftmostRefPos <= _stack[0]->_lastIndex;
}

ColumnIterator::VisitCache *DefaultColumnIterator::getVisitCache() const
{
    return &_visitCache;
}

void DefaultColumnIterator::clearVisitCache() const
{
    for (VisitCache::iterator i = _visitCache.begin();
         i != _visitCache.end(); ++i)
    {
        delete i->second;
    }
    _visitCache.clear();
}

void DefaultColumnIterator::setVisitCache(ColumnIterator::VisitCache *visitCache) const
{
    clearVisitCache();
    _visitCache = *visitCache;
}

void DefaultColumnIterator::print(ostream& os) const
{
  const ColumnIterator::ColumnMap* cmap = getColumnMap();
  for (ColumnIterator::ColumnMap::const_iterator i = cmap->begin();
       i != cmap->end(); ++i)
  {
    os << i->first->getName() << ": ";
    for (size_t j = 0; j < i->second->size(); ++j)
    {
      os << i->second->at(j)->getArrayIndex() << ", ";
    }
    os << "\n";
  }
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
  clearTree();
  _break = false;
  _leftmostRefPos = _stack[0]->_index;

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
      if (_reversed == true)
      {
        topIt->_it->toReverseInPlace();
        topIt->_dna->toReverse();
      }
    }
    // otherwise, we scan forward from last visisted column
    else
    {
      assert(topIt->_it.get() != NULL);
      bool rev = topIt->_it->getReversed();
      assert(rev == _reversed);
      if (rev == true)
      {
        topIt->_it->toReverseInPlace();
      }
      assert(topIt->_it->getReversed() == false);

      // catch up to nextfreeindex
      topIt->_it->slice(0, 0);
      while (topIt->_it->overlaps(_stack.top()->_index) == false)
      {
        topIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - 
                                          topIt->_it->getStartPosition());
      topIt->_it->slice(offset, topIt->_it->getLength() - offset - 1);
      topIt->_dna->jumpTo(_stack.top()->_index);
      if (rev == true)
      {
        assert(topIt->_dna->getReversed() == true);
        topIt->_it->toReverseInPlace();
      }
    }
    assert(topIt->_it->getReversed() == _reversed &&
           topIt->_dna->getReversed() == _reversed);
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
    if (!_onlyOrthologs) {
        updateNextTopDup(topIt);
    }
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
      if (_reversed == true)
      {
        bottomIt->_it->toReverseInPlace();
        bottomIt->_dna->toReverse();
      }
    }
    else
    {
      assert(bottomIt->_it.get() != NULL);
      bool rev = bottomIt->_it->getReversed();
      assert(rev == _reversed);
      if (rev == true)
      {
        bottomIt->_it->toReverseInPlace();
      }
      assert(bottomIt->_it->getReversed() == false);

      // catch up to nextfreeindex
      bottomIt->_it->slice(0, 0);
      while (bottomIt->_it->overlaps(_stack.top()->_index) == false)
      {
        bottomIt->_it->toRight();
      }
      hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - 
                                          bottomIt->_it->getStartPosition());
      bottomIt->_it->slice(offset, bottomIt->_it->getLength() - offset - 1);
      bottomIt->_dna->jumpTo(_stack.top()->_index);
      if (rev == true)
      {
        assert(bottomIt->_dna->getReversed() == true);
        bottomIt->_it->toReverseInPlace();
      }
    }

    assert(bottomIt->_it->getReversed() == _reversed &&
           bottomIt->_dna->getReversed() == _reversed);
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
    if (_top->getEndOffset() == 0 && _top->isLast() == false)
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

// Builds a "gene"-tree node and labels it properly.
static stTree *getTreeNode(SegmentIteratorConstPtr segIt)
{
  // Make sure the segment is sliced to only 1 base.
  assert(segIt->getStartPosition() == segIt->getEndPosition());
  stTree *ret = stTree_construct();
  const Genome *genome = segIt->getGenome();
  const Sequence *seq = genome->getSequenceBySite(segIt->getStartPosition());

  stringstream ss;
  ss << segIt->getGenome()->getName() << "." << seq->getName() << "|" << segIt->getStartPosition() - seq->getStartPosition();
  stTree_setLabel(ret, ss.str().c_str());

  DNAIteratorConstPtr *dnaIt = new DNAIteratorConstPtr(genome->getDNAIterator(segIt->getStartPosition()));
  if (segIt->getReversed()) {
    (*dnaIt)->toReverse();
  }
  stTree_setClientData(ret, (void *) dnaIt);

  return ret;
}

// Recursive part of buildTree
// tree parameter represents node corresponding to the genome with
// bottom segment botIt
static void buildTreeR(BottomSegmentIteratorConstPtr botIt, stTree *tree)
{
  const Genome *genome = botIt->getGenome();

  // attach a node and recurse for each of this segment's children
  // (and paralogous segments)
  for (hal_size_t i = 0; i < botIt->getNumChildren(); i++) {
    if (botIt->hasChild(i)) {
      const Genome *child = genome->getChild(i);
      TopSegmentIteratorConstPtr topIt = child->getTopSegmentIterator();
      topIt->toChild(botIt, i);
      stTree *canonicalParalog = getTreeNode(topIt);
      stTree_setParent(canonicalParalog, tree);
      if (topIt->hasParseDown()) {
        BottomSegmentIteratorConstPtr childBotIt = child->getBottomSegmentIterator();
        childBotIt->toParseDown(topIt);
        buildTreeR(childBotIt, canonicalParalog);
      }
      // Traverse the paralogous segments cycle and add those segments as well
      if (topIt->hasNextParalogy()) {
        topIt->toNextParalogy();
        while(!topIt->isCanonicalParalog()) {
          stTree *paralog = getTreeNode(topIt);
          stTree_setParent(paralog, tree);
          if(topIt->hasParseDown()) {
            BottomSegmentIteratorConstPtr childBotIt = child->getBottomSegmentIterator();
            childBotIt->toParseDown(topIt);
            buildTreeR(childBotIt, paralog);
          }
          topIt->toNextParalogy();
        }
      }
    }
  }
}

// Build a gene-tree from a column iterator.
stTree *DefaultColumnIterator::getTree() const
{
  if (_onlyOrthologs) {
    // Because the tree-finding code goes all the way up the column
    // tree and I'm too lazy to make it smarter.
    throw hal_exception("Cannot get the tree for a column iterator "
                        "which only displays orthologs.");
  }
  if (_tree != NULL) {
    return _tree;
  } else {
    // Get any base from the column to begin building the tree
    const ColumnIterator::ColumnMap *colMap = getColumnMap();
    ColumnIterator::ColumnMap::const_iterator colMapIt = colMap->begin();
    const Sequence *sequence = NULL;
    hal_index_t index = NULL_INDEX;
    while (colMapIt != colMap->end()) {
      if (!colMapIt->second->empty()) {
        // Found a non-empty column map entry, just take the index and
        // sequence of the first base found
        sequence = colMapIt->first;
        index = colMapIt->second->at(0)->getArrayIndex();
        break;
      }
      colMapIt++;
    }
    assert(sequence != NULL && index != NULL_INDEX);
    const Genome *genome = sequence->getGenome();

    // Get the bottom segment that is the common ancestor of all entries
    TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorConstPtr botIt;
    if (genome->getNumTopSegments() == 0) {
      // The reference is the root genome.
      botIt = genome->getBottomSegmentIterator();
      botIt->toSite(index);
    } else {
      // Keep heading up the tree until we hit the root segment.
      topIt->toSite(index);
      while (topIt->hasParent()) {
        const Genome *parent = topIt->getGenome()->getParent();
        botIt = parent->getBottomSegmentIterator();
        botIt->toParent(topIt);
        if(parent->getParent() == NULL || !botIt->hasParseUp()) {
          // Reached root genome
          break;
        }
        topIt = parent->getTopSegmentIterator();
        topIt->toParseUp(botIt);
      }
    }

    stTree *tree = NULL;
    if(topIt->hasParent() == false && topIt->getGenome() == genome && genome->getNumBottomSegments() == 0) {
      // Handle insertions in leaves. botIt doesn't point anywhere since
      // there are no bottom segments.
      tree = getTreeNode(topIt);
    } else {
      tree = getTreeNode(botIt);
      buildTreeR(botIt, tree);
    }

    if (_onlyOrthologs || _noDupes || !_targets.empty()) {
      // The gene tree, at this point, always represents the full
      // induced tree found in the HAL graph. If we are showing part
      // of the full column, we should make sure to give only the
      // corresponding part of the tree.
      //getInducedTree(tree);
    }

    assert(tree != NULL);
    _tree = tree;
    return tree;
  }
}

static void clearTree_R(stTree *tree)
{
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++)
  {
    clearTree_R(stTree_getChild(tree, i));
  }
  delete (DNAIteratorConstPtr *) stTree_getClientData(tree);
}

void DefaultColumnIterator::clearTree() const
{
  if (_tree != NULL)
  {
    clearTree_R(_tree);
    stTree_destruct(_tree);
    _tree = NULL;
  }
}

void DefaultColumnIterator::updateParent(LinkedTopIterator* topIt) const
{
  const Genome* genome = topIt->_it->getTopSegment()->getGenome();
  if (!_break && topIt->_it->hasParent() && parentInScope(genome) &&
      (!_noDupes || topIt->_it->isCanonicalParalog()))
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
  const Genome* genome = bottomIt->_it->getBottomSegment()->getGenome();
  if (!_break && bottomIt->_it->hasChild(index) && childInScope(genome, index))
  {
    assert(index < bottomIt->_children.size());
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
  const Genome* genome =  topIt->_it->getTopSegment()->getGenome();
  if (_break || _noDupes == true ||
      topIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX ||
      genome->getParent() == NULL || parentInScope(genome) == false)
  {
    return;
  }

  hal_index_t firstIndex = topIt->_it->getTopSegment()->getArrayIndex();
  LinkedTopIterator* currentTopIt = topIt;

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
  if (!_break && bottomIt->_it->hasParseUp())
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
    if (!_onlyOrthologs) {
        updateNextTopDup(bottomIt->_topParse);
    }
  }
}
 
void DefaultColumnIterator::updateParseDown(LinkedTopIterator* topIt) const
{
  if (!_break && topIt->_it->hasParseDown())
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

  if (_unique == true || _stack.size() > 1)
  {
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
  }
  _stack.top()->_index = index;
}

bool DefaultColumnIterator::colMapInsert(DNAIteratorConstPtr dnaIt) const
{
  const Sequence* sequence = dnaIt->getSequence();
  const Genome* genome = dnaIt->getGenome();
  assert(sequence != NULL);
  
  // All reference bases need to get added to the cache
  bool updateCache = genome == _stack[0]->_sequence->getGenome();
  if (_maxInsertionLength == 0)
  {
    // Unless we don't do indels.  Here we just add reference elements
    // that are to right of the starting point
    assert (_stack.size() == 1);
    updateCache = 
       updateCache && _stack.top()->_firstIndex < dnaIt->getArrayIndex();
  }
  for (size_t i = 1; i < _stack.size() && !updateCache; ++i)
  {
    if (genome == _stack[i]->_sequence->getGenome())
    {
      updateCache = true;
    }
  }
  // try to avoid building cache if we don't want or need it
  if (_unique == false && _maxInsertionLength == 0)
  {
    updateCache = false;
  }

  bool found = false;
  VisitCache::iterator cacheIt = _visitCache.find(genome);
  if (updateCache == true)
  {
    if (cacheIt == _visitCache.end())
    {
      PositionCache* newSet = new PositionCache();
      cacheIt = _visitCache.insert(pair<const Genome*, PositionCache*>(
                                     genome, newSet)).first;
    }
    found = cacheIt->second->insert(dnaIt->getArrayIndex()) == false;
  }
  else
  {
    found = cacheIt != _visitCache.end() && 
       cacheIt->second->find(dnaIt->getArrayIndex()) == true;
  }

  // insert into the column data structure to pass out to client
  if (found == false && 
      (!_noAncestors || genome->getNumChildren() == 0) &&
      (_targets.empty() || _targets.find(genome) != _targets.end()))
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

  // update leftmost ref pos which is used by isCanonicalOnRef() 
  if (genome == _stack[0]->_sequence->getGenome())
  {
    _leftmostRefPos = min(_leftmostRefPos, dnaIt->getArrayIndex());
  }
  return !found;
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
