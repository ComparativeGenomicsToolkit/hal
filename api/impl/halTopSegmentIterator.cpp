/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halTopSegmentIterator.h"
#include "halBottomSegmentIterator.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE OVERRIDE
//////////////////////////////////////////////////////////////////////////////
std::ostream &TopSegmentIterator::print(ostream &os) const {
    os << "TopSegIt: ";
    SegmentIterator::print(os);

    hal_index_t ai = getArrayIndex();
    bool offRight =
        isTop() ? ai >= (hal_index_t)getGenome()->getNumTopSegments() : ai >= (hal_index_t)getGenome()->getNumBottomSegments();

    if (ai != NULL_INDEX && !offRight) {
        os << " pIdx=" << tseg()->getParentIndex() << " npIdx=" << tseg()->getNextParalogyIndex();
    }
    return os;
}

//////////////////////////////////////////////////////////////////////////////
// TOP SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void TopSegmentIterator::toChild(const BottomSegmentIteratorPtr &botSegIt, hal_size_t child) {
    _topSegment->setArrayIndex(botSegIt->getGenome()->getChild(child), botSegIt->bseg()->getChildIndex(child));
    _startOffset = botSegIt->getStartOffset();
    _endOffset = botSegIt->getEndOffset();
    _reversed = botSegIt->getReversed();
    if (botSegIt->bseg()->getChildReversed(child) == true) {
        toReverse();
    }
    assert(inRange() == true);
}

void TopSegmentIterator::toChildG(const BottomSegmentIteratorPtr &botSegIt, const Genome *childGenome) {
    hal_index_t childIndex = botSegIt->getGenome()->getChildIndex(childGenome);
    assert(childIndex != NULL_INDEX);
    toChild(botSegIt, childIndex);
    assert(inRange() == true);
    assert(getGenome() == childGenome);
}

void TopSegmentIterator::toParseUp(const BottomSegmentIteratorPtr &botSegIt) {
    Genome *genome = botSegIt->getGenome();
    hal_index_t index = botSegIt->bseg()->getTopParseIndex();

    _topSegment->setArrayIndex(genome, index);
    _reversed = botSegIt->getReversed();

    hal_index_t startPos = botSegIt->getStartPosition();

    while (startPos >= _topSegment->getStartPosition() + (hal_index_t)_topSegment->getLength()) {
        _topSegment->setArrayIndex(genome, ++index);
    }

    if (_reversed == false) {
        _startOffset = startPos - _topSegment->getStartPosition();
        hal_index_t topEnd = _topSegment->getStartPosition() + (hal_index_t)_topSegment->getLength();
        hal_index_t botEnd = botSegIt->getStartPosition() + (hal_index_t)botSegIt->getLength();
        _endOffset = max((hal_index_t)0, topEnd - botEnd);
    } else {
        _startOffset = _topSegment->getStartPosition() + _topSegment->getLength() - 1 - startPos;
        hal_index_t topEnd = _topSegment->getStartPosition();
        hal_index_t botEnd = botSegIt->getStartPosition() - (hal_index_t)botSegIt->getLength() + 1;
        _endOffset = max((hal_index_t)0, botEnd - topEnd);
    }
    assert(_startOffset + _endOffset <= _topSegment->getLength());
    assert(inRange() == true);
}

TopSegmentIteratorPtr TopSegmentIterator::clone() const {
    TopSegmentIteratorPtr newIt = getGenome()->getTopSegmentIterator(getArrayIndex());
    if (_reversed) {
        newIt->toReverse();
    }
    newIt->slice(_startOffset, _endOffset);
    return newIt;
}

void TopSegmentIterator::copy(const TopSegmentIteratorPtr &topSegIt) {
    _topSegment->setArrayIndex(topSegIt->getGenome(), topSegIt->getArrayIndex());
    _startOffset = topSegIt->getStartOffset();
    _endOffset = topSegIt->getEndOffset();
    _reversed = topSegIt->getReversed();
}

void TopSegmentIterator::toNextParalogy() {
    assert(_topSegment->getNextParalogyIndex() != NULL_INDEX);
    assert(_topSegment->getNextParalogyIndex() != _topSegment->getArrayIndex());
    bool rev = _topSegment->getParentReversed();
    _topSegment->setArrayIndex(getGenome(), _topSegment->getNextParalogyIndex());
    if (_topSegment->getParentReversed() != rev) {
        toReverse();
    }
}
