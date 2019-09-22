/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halColumnIterator.h"
#include "halBottomSegmentIterator.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

ColumnIterator::ColumnIterator(const Genome *reference, const set<const Genome *> *targets, hal_index_t columnIndex,
                               hal_index_t lastColumnIndex, hal_size_t maxInsertLength, bool noDupes, bool noAncestors,
                               bool reverseStrand, bool unique, bool onlyOrthologs)
    : _maxInsertionLength(maxInsertLength), _noDupes(noDupes), _noAncestors(noAncestors),
      _treeCache(NULL), _unique(unique), _onlyOrthologs(onlyOrthologs) {
    assert(columnIndex >= 0 && lastColumnIndex >= columnIndex && lastColumnIndex < (hal_index_t)reference->getSequenceLength());

    // allocate temp iterators
    if (reference->getNumTopSegments() > 0) {
        _top = reference->getTopSegmentIterator(0);
        _next = _top->clone();
    } else if (reference->getChild(0) != NULL) {
        _top = reference->getChild(0)->getTopSegmentIterator(0);
        _next = _top->clone();
    }

    if (_maxInsertionLength > 0) {
        // need to allocate the rearrangement from
        if (reference->getParent() != NULL) {
            _rearrangement = reference->getRearrangement(0, 0, 1., true);
        } else if (reference->getNumChildren() > 0) {
            _rearrangement = reference->getChild(0)->getRearrangement(0, 0, 1., true);
        }
    }
    const Sequence *sequence = reference->getSequenceBySite(columnIndex);
    assert(sequence != NULL);
    _ref = sequence;

    // compute all genomes in search scope (spanning tree of reference and targets)
    // if targets is empty we just visit everything.
    if (targets != NULL && !targets->empty()) {
        _targets = *targets;
        _targets.insert(reference);
        getGenomesInSpanningTree(_targets, _scope);
    }

    // note columnIndex in genome (not sequence) coordinates
    _stack.push(sequence, columnIndex, lastColumnIndex, reverseStrand);

    toRight();
}

ColumnIterator::~ColumnIterator() {
    eraseColMap();
    clearVisitCache();
    clearTree();
}

void ColumnIterator::toRight() {
    clearTree();

    // keep the current position so that when client calls
    // getReferenceXXX() methods, they get the state before
    // toRight is called.
    _prevRefSequence = _ref;
    _prevRefIndex = _stack[0]->_index - _ref->getStartPosition();

    // compatible with old interface which allowed toRight() to go out
    // of bounds without crashing.
    if (_stack.size() == 1 && !_stack.topInBounds()) {
        return;
    }
    assert(_indelStack.size() == 0);

    do {
        // clean stack
        nextFreeIndex();
        while (_stack.size() > 1 && !_stack.topInBounds()) {
            _stack.popDelete();
            nextFreeIndex();
        }

        // compatible with old interface which allowed toRight() to go out
        // of bounds without crashing.
        if (_stack.size() == 1 && !_stack.topInBounds()) {
            return;
        }

        _indelStack.clear();

        bool init = _stack.top()->_index == _stack.top()->_firstIndex ||
                    (_stack.top()->_bottom._it.get() == NULL && _stack.top()->_top._it.get() == NULL);

        recursiveUpdate(init);

        // move the index right
        if (_stack.top()->_reversed) {
            _stack.top()->_index--;
        } else {
            _stack.top()->_index++;
        }

        // jump to next sequence in genome if necessary
        const Sequence *seq = _stack.top()->_sequence;
        if (_stack.size() == 1 &&
            (_stack.top()->_index < seq->getStartPosition() ||
             (_stack.top()->_index >= (hal_index_t)(seq->getStartPosition() + seq->getSequenceLength()) &&
            _stack.top()->_index < (hal_index_t)(seq->getGenome()->getSequenceLength())))) {
            _stack.top()->_sequence = seq->getGenome()->getSequenceBySite(_stack.top()->_index);
            assert(_stack.top()->_sequence != NULL);
            _ref = _stack.top()->_sequence;
        }
    } while (_break == true);

    // push the indel stack.
    _stack.pushStack(_indelStack);

    // clean stack again
    nextFreeIndex();
    while (_stack.size() > 1 && !_stack.topInBounds()) {
        _stack.popDelete();
        nextFreeIndex();
    }

#ifndef NDEBUG
    set<pair<const Sequence *, hal_index_t>> coordSet;
    ColumnMap::const_iterator i, iNext;
    DNASet::const_iterator j;
    for (i = _colMap.begin(); i != _colMap.end(); ++i) {
        // check that the same coordinate not present for the same sequence
        for (j = i->second->begin(); j != i->second->end(); ++j) {
            pair<const Sequence *, hal_index_t> data(i->first, (*j)->getArrayIndex());
            assert(coordSet.insert(data).second == true);
        }
    }
#endif
}

void ColumnIterator::toSite(hal_index_t columnIndex, hal_index_t lastColumnIndex, bool clearCache) {
    clearTree();

    const Genome *reference = getReferenceGenome();
    assert(columnIndex >= 0 && lastColumnIndex >= columnIndex && lastColumnIndex < (hal_index_t)reference->getSequenceLength());

    const Sequence *sequence = reference->getSequenceBySite(columnIndex);
    assert(sequence != NULL);
    _ref = sequence;
    _stack.clear();
    _indelStack.clear();
    if (clearCache == true) {
        clearVisitCache();
    }
    defragment();
    // note columnIndex in genome (not sequence) coordinates
    _stack.push(sequence, columnIndex, lastColumnIndex);
    toRight();
    assert(getReferenceSequencePosition() + sequence->getStartPosition() == columnIndex);
}

bool ColumnIterator::lastColumn() const {
    return _stack.size() == 1 && _stack.top()->pastEnd();
}

const Genome *ColumnIterator::getReferenceGenome() const {
    return _prevRefSequence->getGenome();
}

const Sequence *ColumnIterator::getReferenceSequence() const {
    return _prevRefSequence;
}

hal_index_t ColumnIterator::getReferenceSequencePosition() const {
    return _prevRefIndex;
}

const ColumnIterator::ColumnMap *ColumnIterator::getColumnMap() const {
    return &_colMap;
}

hal_index_t ColumnIterator::getArrayIndex() const {
    assert(_stack.size() > 0);
    return _stack[0]->_index;
}

void ColumnIterator::defragment() {
    ColumnMap::iterator i = _colMap.begin();
    ColumnMap::iterator next;
    while (i != _colMap.end()) {
        next = i;
        ++next;
        if (i->second->empty()) {
            delete i->second;
            _colMap.erase(i);
        }
        i = next;
    }

    _stack.resetLinks();
}

bool ColumnIterator::isCanonicalOnRef() const {
    assert(_stack.size() > 0);
    assert(_leftmostRefPos >= 0 && (hal_size_t)_leftmostRefPos < _stack[0]->_sequence->getGenome()->getSequenceLength());
    return _leftmostRefPos >= _stack[0]->_firstIndex && _leftmostRefPos <= _stack[0]->_lastIndex;
}

ColumnIterator::VisitCache *ColumnIterator::getVisitCache() {
    return &_visitCache;
}

void ColumnIterator::clearVisitCache() {
    for (VisitCache::iterator i = _visitCache.begin(); i != _visitCache.end(); ++i) {
        delete i->second;
    }
    _visitCache.clear();
}

void ColumnIterator::setVisitCache(ColumnIterator::VisitCache *visitCache) {
    clearVisitCache();
    _visitCache = *visitCache;
}

void ColumnIterator::print(ostream &os) const {
    const ColumnIterator::ColumnMap *cmap = getColumnMap();
    for (ColumnIterator::ColumnMap::const_iterator i = cmap->begin(); i != cmap->end(); ++i) {
        os << i->first->getName() << ": ";
        for (size_t j = 0; j < i->second->size(); ++j) {
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
void ColumnIterator::recursiveUpdate(bool init) {
    resetColMap();
    clearTree();
    _break = false;
    _leftmostRefPos = _stack[0]->_index;

    const Sequence *refSequence = _stack.top()->_sequence;
    const Genome *refGenome = refSequence->getGenome();
    if (refSequence->getNumTopSegments() > 0) {
        assert(_stack.size() > 0);
        LinkedTopIterator *linkTopIt = &_stack.top()->_top;
        // first column, we search the genome for the site
        if (init == true) {
            linkTopIt->_it = refSequence->getTopSegmentIterator();
            linkTopIt->_it->toSite(_stack.top()->_index, true);
            linkTopIt->_dna = refGenome->getDnaIterator(_stack.top()->_index);
            if (_stack.top()->_reversed) {
                linkTopIt->_it->toReverseInPlace();
                linkTopIt->_dna->toReverse();
            }
        }
        // otherwise, we scan forward from last visisted column
        else {
            assert(linkTopIt->_it.get() != NULL);

            // catch up to nextfreeindex
            linkTopIt->_it->slice(0, 0);
            while (linkTopIt->_it->overlaps(_stack.top()->_index) == false) {
                linkTopIt->_it->toRight();
            }
            bool rev = linkTopIt->_it->getReversed();
            if (rev == true) {
                linkTopIt->_it->toReverseInPlace();
            }
            assert(linkTopIt->_it->getReversed() == false);
            hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - linkTopIt->_it->getStartPosition());
            linkTopIt->_it->slice(offset, linkTopIt->_it->getLength() - offset - 1);
            linkTopIt->_dna->jumpTo(_stack.top()->_index);
            if (rev == true) {
                assert(linkTopIt->_dna->getReversed() == true);
                linkTopIt->_it->toReverseInPlace();
            }
        }
        assert(linkTopIt->_it->getStartPosition() == linkTopIt->_dna->getArrayIndex());
        assert(linkTopIt->_dna->getArrayIndex() == _stack.top()->_index);
        assert(_stack.top()->_index <= _stack.top()->_lastIndex);
        assert(linkTopIt->_it->getStartPosition() == linkTopIt->_dna->getArrayIndex());

        if (colMapInsert(linkTopIt->_dna) == false) {
            _break = true;
            return;
        }
        handleDeletion(linkTopIt->_it);
        updateParent(linkTopIt);
        if (!_onlyOrthologs) {
            updateNextTopDup(linkTopIt);
        }
        updateParseDown(linkTopIt);
    }

    else {
        assert(_stack.size() > 0);
        LinkedBottomIterator *linkBotIt = &_stack.top()->_bottom;
        if (init == true) {
            linkBotIt->_it = refSequence->getBottomSegmentIterator();
            linkBotIt->_it->toSite(_stack.top()->_index, true);
            linkBotIt->_dna = refGenome->getDnaIterator(_stack.top()->_index);
            if (_stack.top()->_reversed) {
                linkBotIt->_it->toReverseInPlace();
                linkBotIt->_dna->toReverse();
            }
        } else {
            assert(linkBotIt->_it.get() != NULL);

            // catch up to nextfreeindex
            linkBotIt->_it->slice(0, 0);
            while (linkBotIt->_it->overlaps(_stack.top()->_index) == false) {
                linkBotIt->_it->toRight();
            }
            bool rev = linkBotIt->_it->getReversed();
            if (rev == true) {
                linkBotIt->_it->toReverseInPlace();
            }
            assert(linkBotIt->_it->getReversed() == false);
            hal_size_t offset = (hal_size_t)abs(_stack.top()->_index - linkBotIt->_it->getStartPosition());
            linkBotIt->_it->slice(offset, linkBotIt->_it->getLength() - offset - 1);
            linkBotIt->_dna->jumpTo(_stack.top()->_index);
            if (rev == true) {
                assert(linkBotIt->_dna->getReversed() == true);
                linkBotIt->_it->toReverseInPlace();
            }
        }

        assert(linkBotIt->_it->getStartPosition() == linkBotIt->_dna->getArrayIndex());
        assert(linkBotIt->_dna->getArrayIndex() == _stack.top()->_index);

        if (colMapInsert(linkBotIt->_dna) == false) {
            _break = true;
            return;
        }
        hal_size_t numChildren = refSequence->getGenome()->getNumChildren();
        if (numChildren > linkBotIt->_children.size()) {
            linkBotIt->_children.resize(numChildren, NULL);
        }
        assert(linkBotIt->_it->getStartPosition() == linkBotIt->_dna->getArrayIndex());
        for (size_t child = 0; child < numChildren; ++child) {
            updateChild(linkBotIt, child);
        }
    }
}

bool ColumnIterator::handleDeletion(const TopSegmentIteratorPtr &inputTopSegIt) {
    if (_maxInsertionLength > 0 && inputTopSegIt->tseg()->hasParent() == true) {
        _top->copy(inputTopSegIt);
        bool reversed = _top->getReversed();
        if (reversed == true) {
            // We need to identify the deletion from the *left* breakpoint on
            // the positive strand. That means the *right* breakpoint when on
            // the negative strand. Therefore, we go right (i.e. left according
            // to the positive strand) when checking the breakpoint.
            _top->toRight();
            _top->toReverse();
        }
        // only handle a deletion if we are immediately left of the breakpoint
        if (_top->getEndOffset() == 0) {
            const Genome *genome = _top->getTopSegment()->getGenome();
            const Genome *parent = genome->getParent();
            _top->slice(0, 0);
            assert(_rearrangement->getAtomic() == true);
            if (_rearrangement->identifyDeletionFromLeftBreakpoint(_top) == true &&
                _rearrangement->getLength() + _stack.top()->_cumulativeSize <= _maxInsertionLength) {
                pair<hal_index_t, hal_index_t> deletedRange = _rearrangement->getDeletedRange();
                assert((hal_size_t)(deletedRange.second - deletedRange.first) == _rearrangement->getLength() - 1);

                BottomSegmentIteratorPtr botSegIt = parent->getBottomSegmentIterator(0);
                botSegIt->toParent(inputTopSegIt);
                _indelStack.push(botSegIt->getBottomSegment()->getSequence(), deletedRange.first, deletedRange.second, botSegIt->getReversed());

                return true;
            }
        }
    }
    return false;
}

bool ColumnIterator::handleInsertion(const TopSegmentIteratorPtr &inputTopSegIt) {
    if (_maxInsertionLength > 0 && inputTopSegIt->tseg()->hasParent() == true) {
        _top->copy(inputTopSegIt);
        bool reversed = _top->getReversed();
        // only handle an insertion if we are immediately left of the break
        if (_top->getEndOffset() == 0 && _top->isLast() == false) {
            _rearrangement->setAtomic(true);
            _top->slice(0, 0);
            _top->toRight();
            if (reversed == true) {
                _top->toReverse();
            }
            assert(_rearrangement->getAtomic() == true);
            if (_rearrangement->identifyInsertionFromLeftBreakpoint(_top) == true &&
                _rearrangement->getLength() + _stack.top()->_cumulativeSize <= _maxInsertionLength) {
                pair<hal_index_t, hal_index_t> insertedRange = _rearrangement->getInsertedRange();
                assert((hal_size_t)(insertedRange.second - insertedRange.first) == _rearrangement->getLength() - 1);

                _indelStack.push(_top->getTopSegment()->getSequence(), insertedRange.first, insertedRange.second, reversed);
            }
        }
    }
    return false;
}

// Builds a "gene"-tree node and labels it properly.
static stTree *getTreeNode(SegmentIteratorPtr segIt) {
    // Make sure the segment is sliced to only 1 base.
    assert(segIt->getStartPosition() == segIt->getEndPosition());
    stTree *ret = stTree_construct();
    const Genome *genome = segIt->getGenome();
    const Sequence *seq = genome->getSequenceBySite(segIt->getStartPosition());

    std::string label(segIt->getGenome()->getName() + "." + seq->getName() + "|" +
                      std::to_string(segIt->getStartPosition() - seq->getStartPosition()));
    stTree_setLabel(ret, label.c_str());

    DnaIteratorPtr *dnaIt = new DnaIteratorPtr(genome->getDnaIterator(segIt->getStartPosition()));
    if (segIt->getReversed()) {
        (*dnaIt)->toReverse();
    }
    stTree_setClientData(ret, (void *)dnaIt);

    return ret;
}

// Recursive part of buildTree
// tree parameter represents node corresponding to the genome with
// bottom segment botIt
static void buildTreeR(BottomSegmentIteratorPtr botSegIt, stTree *tree) {
    const Genome *genome = botSegIt->getGenome();

    // attach a node and recurse for each of this segment's children
    // (and paralogous segments)
    for (hal_size_t i = 0; i < botSegIt->bseg()->getNumChildren(); i++) {
        if (botSegIt->bseg()->hasChild(i)) {
            const Genome *child = genome->getChild(i);
            TopSegmentIteratorPtr topSegIt = child->getTopSegmentIterator();
            topSegIt->toChild(botSegIt, i);
            stTree *canonicalParalog = getTreeNode(topSegIt);
            stTree_setParent(canonicalParalog, tree);
            if (topSegIt->tseg()->hasParseDown()) {
                BottomSegmentIteratorPtr childBotSegIt = child->getBottomSegmentIterator();
                childBotSegIt->toParseDown(topSegIt);
                buildTreeR(childBotSegIt, canonicalParalog);
            }
            // Traverse the paralogous segments cycle and add those segments as well
            if (topSegIt->tseg()->hasNextParalogy()) {
                topSegIt->toNextParalogy();
                while (!topSegIt->tseg()->isCanonicalParalog()) {
                    stTree *paralog = getTreeNode(topSegIt);
                    stTree_setParent(paralog, tree);
                    if (topSegIt->tseg()->hasParseDown()) {
                        BottomSegmentIteratorPtr childBotSegIt = child->getBottomSegmentIterator();
                        childBotSegIt->toParseDown(topSegIt);
                        buildTreeR(childBotSegIt, paralog);
                    }
                    topSegIt->toNextParalogy();
                }
            }
        }
    }
}

// Build cached gene-tree from a column iterator.
stTree *ColumnIterator::buildTree() const {
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
    TopSegmentIteratorPtr topSegIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botSegIt;
    if (genome->getNumTopSegments() == 0) {
        // The reference is the root genome.
        botSegIt = genome->getBottomSegmentIterator();
        botSegIt->toSite(index);
    } else {
        // Keep heading up the tree until we hit the root segment.
        topSegIt->toSite(index);
        while (topSegIt->tseg()->hasParent()) {
            const Genome *parent = topSegIt->getGenome()->getParent();
            botSegIt = parent->getBottomSegmentIterator();
            botSegIt->toParent(topSegIt);
            if (parent->getParent() == NULL || !botSegIt->bseg()->hasParseUp()) {
                // Reached root genome
                break;
            }
            topSegIt = parent->getTopSegmentIterator();
            topSegIt->toParseUp(botSegIt);
        }
    }

    stTree *tree = NULL;
    if ((genome->getNumTopSegments() != 0) && (not topSegIt->tseg()->hasParent()) && (topSegIt->getGenome() == genome) &&
        (genome->getNumBottomSegments() == 0)) {
        // Handle insertions in leaves. botSegIt doesn't point anywhere since
        // there are no bottom segments.
        tree = getTreeNode(topSegIt);
    } else {
        tree = getTreeNode(botSegIt);
        buildTreeR(botSegIt, tree);
    }

    if (_onlyOrthologs || _noDupes || !_targets.empty()) {
        // The gene tree, at this point, always represents the full
        // induced tree found in the HAL graph. If we are showing part
        // of the full column, we should make sure to give only the
        // corresponding part of the tree.
        // getInducedTree(tree);
    }

    assert(tree != NULL);
    return tree;
}

// Build a gene-tree from a column iterator.
stTree *ColumnIterator::getTree() const {
    if (_onlyOrthologs || _noDupes) {
        // Because the tree-finding code goes all the way up the column
        // tree and I'm too lazy to make it smarter.
        throw hal_exception("Cannot get the tree for a column iterator "
                            "which only displays orthologs.");
    }
    if (_treeCache == NULL) {
        _treeCache = buildTree();
    }
    return _treeCache;
}

static void clearTree_R(stTree *tree) {
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        clearTree_R(stTree_getChild(tree, i));
    }
    delete (DnaIteratorPtr *)stTree_getClientData(tree);
}

void ColumnIterator::clearTree() {
    if (_treeCache != NULL) {
        clearTree_R(_treeCache);
        stTree_destruct(_treeCache);
        _treeCache = NULL;
    }
}

void ColumnIterator::updateParent(LinkedTopIterator *linkTopIt) {
    const Genome *genome = linkTopIt->_it->getTopSegment()->getGenome();
    if (!_break && linkTopIt->_it->tseg()->hasParent() && parentInScope(genome) &&
        (!_noDupes || linkTopIt->_it->tseg()->isCanonicalParalog())) {
        const Genome *parentGenome = genome->getParent();

        // no linked iterator for parent.  we create a new one and add the
        // link in both directions
        if (linkTopIt->_parent == NULL) {
            assert(parentGenome != NULL);
            linkTopIt->_parent = linkTopIt->_entry->newBottom();
            linkTopIt->_parent->_it = parentGenome->getBottomSegmentIterator();
            linkTopIt->_parent->_dna = parentGenome->getDnaIterator();
            hal_size_t numChildren = parentGenome->getNumChildren();
            if (numChildren > linkTopIt->_parent->_children.size()) {
                linkTopIt->_parent->_children.resize(numChildren, NULL);
            }
            for (hal_size_t i = 0; i < numChildren; ++i) {
                if (parentGenome->getChild(i) == genome) {
                    linkTopIt->_parent->_children[i] = linkTopIt;
                }
            }
        }

        // advance the parent's iterator to match linkTopIt's (which should
        // already have been updated.
        linkTopIt->_parent->_it->toParent(linkTopIt->_it);
        linkTopIt->_parent->_dna->jumpTo(linkTopIt->_parent->_it->getStartPosition());
        linkTopIt->_parent->_dna->setReversed(linkTopIt->_parent->_it->getReversed());
        if (colMapInsert(linkTopIt->_parent->_dna) == false) {
            _break = true;
            return;
        }

        // recurse on parent's parse edge
        updateParseUp(linkTopIt->_parent);
        if (linkTopIt->_parent->_it->bseg()->hasParseUp() && linkTopIt->_parent->_topParse->_it.get() != NULL) {
            handleDeletion(linkTopIt->_parent->_topParse->_it);
        }

        // recurse on parent's child edges (siblings to linkTopIt)
        for (hal_size_t i = 0; i < linkTopIt->_parent->_children.size(); ++i) {
            if (linkTopIt->_parent->_children[i] == NULL ||
                linkTopIt->_parent->_children[i]->_it->getTopSegment()->getGenome() != genome) {
                updateChild(linkTopIt->_parent, i);
            }
        }
    }
}

void ColumnIterator::updateChild(LinkedBottomIterator *linkBotIt, hal_size_t index) {
    const Genome *genome = linkBotIt->_it->getBottomSegment()->getGenome();
    if (!_break && linkBotIt->_it->bseg()->hasChild(index) && childInScope(genome, index)) {
        assert(index < linkBotIt->_children.size());
        const Genome *childGenome = genome->getChild(index);

        // no linked iterator for child. we create a new one and add linke in
        // both directions
        if (linkBotIt->_children[index] == NULL) {
            assert(childGenome != NULL);
            linkBotIt->_children[index] = linkBotIt->_entry->newTop();
            linkBotIt->_children[index]->_it = childGenome->getTopSegmentIterator();
            linkBotIt->_children[index]->_dna = childGenome->getDnaIterator();
            linkBotIt->_children[index]->_parent = linkBotIt;
        }

        // advance the child's iterator to match linkBotIt's (which should
        // have already been updated)
        linkBotIt->_children[index]->_it->toChild(linkBotIt->_it, index);
        linkBotIt->_children[index]->_dna->jumpTo(linkBotIt->_children[index]->_it->getStartPosition());
        linkBotIt->_children[index]->_dna->setReversed(linkBotIt->_children[index]->_it->getReversed());
        if (colMapInsert(linkBotIt->_children[index]->_dna) == false) {
            _break = true;
            return;
        }
        handleInsertion(linkBotIt->_children[index]->_it);

        // recurse on paralgous siblings
        updateNextTopDup(linkBotIt->_children[index]);

        // recurse on child's parse edge
        updateParseDown(linkBotIt->_children[index]);
    }
}

void ColumnIterator::updateNextTopDup(LinkedTopIterator *linkTopIt) {
    assert(linkTopIt->_it.get() != NULL);
    const Genome *genome = linkTopIt->_it->getTopSegment()->getGenome();
    if (_break || _noDupes == true || linkTopIt->_it->getTopSegment()->getNextParalogyIndex() == NULL_INDEX ||
        genome->getParent() == NULL || parentInScope(genome) == false) {
        return;
    }

    hal_index_t firstIndex = linkTopIt->_it->getTopSegment()->getArrayIndex();
    LinkedTopIterator *currentTopIt = linkTopIt;

    do {
        // no linked iterator for paralog. we create a new one and add link
        if (currentTopIt->_nextDup == NULL) {
            currentTopIt->_nextDup = currentTopIt->_entry->newTop();
            currentTopIt->_nextDup->_it = genome->getTopSegmentIterator();
            currentTopIt->_nextDup->_dna = genome->getDnaIterator();
            currentTopIt->_nextDup->_parent = currentTopIt->_parent;
        }

        // advance the dups's iterator to match currentTopIt's (which should
        // have already been updated)
        currentTopIt->_nextDup->_it = currentTopIt->_it->clone();
        currentTopIt->_nextDup->_it->toNextParalogy();
        currentTopIt->_nextDup->_dna->jumpTo(currentTopIt->_nextDup->_it->getStartPosition());
        currentTopIt->_nextDup->_dna->setReversed(currentTopIt->_nextDup->_it->getReversed());
        if (colMapInsert(currentTopIt->_nextDup->_dna) == false) {
            _break = true;
            return;
        }
        handleInsertion(currentTopIt->_nextDup->_it);

        // recurse on duplicate's parse edge
        updateParseDown(currentTopIt->_nextDup);

        // advance current it to the next paralog
        currentTopIt = currentTopIt->_nextDup;
    } while (currentTopIt->_it->getTopSegment()->getNextParalogyIndex() != NULL_INDEX &&
             currentTopIt->_it->getTopSegment()->getNextParalogyIndex() != firstIndex);
}

void ColumnIterator::updateParseUp(LinkedBottomIterator *linkBotIt) {
    if (!_break && linkBotIt->_it->bseg()->hasParseUp()) {
        const Genome *genome = linkBotIt->_it->getBottomSegment()->getGenome();

        // no linked iterator for top parse, we create a new one
        if (linkBotIt->_topParse == NULL) {
            linkBotIt->_topParse = linkBotIt->_entry->newTop();
            linkBotIt->_topParse->_it = genome->getTopSegmentIterator();
            linkBotIt->_topParse->_dna = genome->getDnaIterator();
            linkBotIt->_topParse->_bottomParse = linkBotIt;
        }

        // advance the parse link's iterator to match linkBotIt
        linkBotIt->_topParse->_it->toParseUp(linkBotIt->_it);
        linkBotIt->_topParse->_dna->jumpTo(linkBotIt->_topParse->_it->getStartPosition());
        linkBotIt->_topParse->_dna->setReversed(linkBotIt->_topParse->_it->getReversed());
        assert(linkBotIt->_topParse->_dna->getArrayIndex() == linkBotIt->_dna->getArrayIndex());

        // recurse on parse link's parent
        updateParent(linkBotIt->_topParse);

        // recurse on parse link's paralogous siblings
        if (!_onlyOrthologs) {
            updateNextTopDup(linkBotIt->_topParse);
        }
    }
}

void ColumnIterator::updateParseDown(LinkedTopIterator *linkTopIt) {
    if (!_break && linkTopIt->_it->tseg()->hasParseDown()) {
        const Genome *genome = linkTopIt->_it->getTopSegment()->getGenome();

        // no linked iterator for down parse, we create a new one
        if (linkTopIt->_bottomParse == NULL) {
            linkTopIt->_bottomParse = linkTopIt->_entry->newBottom();
            linkTopIt->_bottomParse->_it = genome->getBottomSegmentIterator();
            linkTopIt->_bottomParse->_dna = genome->getDnaIterator();
            linkTopIt->_bottomParse->_topParse = linkTopIt;
            hal_size_t numChildren = genome->getNumChildren();
            if (numChildren > linkTopIt->_bottomParse->_children.size()) {
                linkTopIt->_bottomParse->_children.resize(numChildren, NULL);
            }
            for (hal_size_t i = 0; i < numChildren; ++i) {
                if (genome->getChild(i) == genome) {
                    linkTopIt->_bottomParse->_children[i] = linkTopIt;
                }
            }
        }

        // advance the parse link's iterator to match linkTopIt
        linkTopIt->_bottomParse->_it->toParseDown(linkTopIt->_it);

        linkTopIt->_bottomParse->_dna->jumpTo(linkTopIt->_bottomParse->_it->getStartPosition());
        linkTopIt->_bottomParse->_dna->setReversed(linkTopIt->_bottomParse->_it->getReversed());
        assert(linkTopIt->_bottomParse->_dna->getArrayIndex() == linkTopIt->_dna->getArrayIndex());

        // recurse on all the link's children
        for (hal_size_t i = 0; i < linkTopIt->_bottomParse->_children.size(); ++i) {
            updateChild(linkTopIt->_bottomParse, i);
        }
    }
}

// moves index "right" until unvisited base is found
// if none exists in the current range, index is left one
// spot out of bounds (invalid) and return false.
void ColumnIterator::nextFreeIndex() {
    hal_index_t index = _stack.top()->_index;

    if (_unique == true || _stack.size() > 1) {
        const VisitCache::iterator cacheIt = _visitCache.find(_stack.top()->_sequence->getGenome());
        if (cacheIt != _visitCache.end()) {
            PositionCache *posCache = cacheIt->second;
            bool found = posCache->find(index);
            while (found == true && index <= _stack.top()->_lastIndex) {
                ++index;
                found = posCache->find(index);
            }
        }
    }
    _stack.top()->_index = index;
}

bool ColumnIterator::colMapInsert(DnaIteratorPtr dnaIt) {
    const Sequence *sequence = dnaIt->getSequence();
    const Genome *genome = dnaIt->getGenome();
    assert(sequence != NULL);

    // All reference bases need to get added to the cache
    bool updateCache = genome == _stack[0]->_sequence->getGenome();
    if (_maxInsertionLength == 0) {
        // Unless we don't do indels.  Here we just add reference elements
        // that are to right of the starting point
        assert(_stack.size() == 1);
        updateCache = updateCache && _stack.top()->_firstIndex < dnaIt->getArrayIndex();
    }
    for (size_t i = 1; i < _stack.size() && !updateCache; ++i) {
        if (genome == _stack[i]->_sequence->getGenome()) {
            updateCache = true;
        }
    }
    // try to avoid building cache if we don't want or need it
    if (_unique == false && _maxInsertionLength == 0) {
        updateCache = false;
    }

    bool found = false;
    VisitCache::iterator cacheIt = _visitCache.find(genome);
    if (updateCache == true) {
        if (cacheIt == _visitCache.end()) {
            PositionCache *newSet = new PositionCache();
            cacheIt = _visitCache.insert(pair<const Genome *, PositionCache *>(genome, newSet)).first;
        }
        found = cacheIt->second->insert(dnaIt->getArrayIndex()) == false;
    } else {
        found = cacheIt != _visitCache.end() && cacheIt->second->find(dnaIt->getArrayIndex()) == true;
    }

    // insert into the column data structure to pass out to client
    if (found == false && (!_noAncestors || genome->getNumChildren() == 0) &&
        (_targets.empty() || _targets.find(genome) != _targets.end())) {
        ColumnMap::iterator i = _colMap.lower_bound(sequence);
        if (i != _colMap.end() && !(_colMap.key_comp()(sequence, i->first))) {
            i->second->push_back(dnaIt);
        } else {
            DNASet *dnaSet = new DNASet();
            dnaSet->push_back(dnaIt);
            _colMap.insert(i, ColumnMap::value_type(sequence, dnaSet));
        }
    }

    // update leftmost ref pos which is used by isCanonicalOnRef()
    if (genome == _stack[0]->_sequence->getGenome()) {
        _leftmostRefPos = min(_leftmostRefPos, dnaIt->getArrayIndex());
    }
    return !found;
}

void ColumnIterator::resetColMap() {
    for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i) {
        i->second->clear();
    }
}

void ColumnIterator::eraseColMap() {
    for (ColumnMap::iterator i = _colMap.begin(); i != _colMap.end(); ++i) {
        delete i->second;
    }
    _colMap.clear();
}
