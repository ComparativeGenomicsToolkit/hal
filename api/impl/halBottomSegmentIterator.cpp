/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halBottomSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

//////////////////////////////////////////////////////////////////////////////
// SEGMENT INTERFACE OVERRIDE
//////////////////////////////////////////////////////////////////////////////
std::ostream &BottomSegmentIterator::print(ostream &os) const {
    os << "BotSegIt: ";
    SegmentIterator::print(os);

    hal_index_t ai = getArrayIndex();
    bool offRight =
        isTop() ? ai >= (hal_index_t)getGenome()->getNumTopSegments() : ai >= (hal_index_t)getGenome()->getNumBottomSegments();

    if (ai != NULL_INDEX && !offRight) {
        os << " numChilds=" << bseg()->getNumChildren();
        for (hal_size_t i = 0; i < bseg()->getNumChildren(); ++i) {
            os << " cI[" << i << "]=" << bseg()->getChildIndex(i);
            os << " cR[" << i << "]=" << bseg()->getChildReversed(i);
        }
    }
    return os;
}

//////////////////////////////////////////////////////////////////////////////
// BOTTOM SEGMENT ITERATOR INTERFACE
//////////////////////////////////////////////////////////////////////////////
void BottomSegmentIterator::toParent(const TopSegmentIteratorPtr &topSegIt) {
    _bottomSegment->setArrayIndex(topSegIt->getGenome()->getParent(), topSegIt->tseg()->getParentIndex());
    _startOffset = topSegIt->getStartOffset();
    _endOffset = topSegIt->getEndOffset();
    _reversed = topSegIt->getReversed();
    if (topSegIt->tseg()->getParentReversed() == true) {
        toReverse();
    }
    assert(inRange() == true);
}

void BottomSegmentIterator::toParseDown(const TopSegmentIteratorPtr &topSegIt) {
    Genome *genome = topSegIt->getGenome();
    hal_index_t index = topSegIt->tseg()->getBottomParseIndex();

    _bottomSegment->setArrayIndex(genome, index);
    _reversed = topSegIt->getReversed();

    hal_index_t startPos = topSegIt->getStartPosition();
    while (startPos >= _bottomSegment->getStartPosition() + (hal_index_t)_bottomSegment->getLength()) {
        _bottomSegment->setArrayIndex(genome, ++index);
    }

    if (_reversed == false) {
        _startOffset = startPos - _bottomSegment->getStartPosition();
        hal_index_t topEnd = _bottomSegment->getStartPosition() + (hal_index_t)_bottomSegment->getLength();
        hal_index_t botEnd = topSegIt->getStartPosition() + (hal_index_t)topSegIt->getLength();
        _endOffset = max((hal_index_t)0, topEnd - botEnd);
    } else {
        _startOffset = _bottomSegment->getStartPosition() + _bottomSegment->getLength() - 1 - startPos;
        hal_index_t topEnd = _bottomSegment->getStartPosition();
        hal_index_t botEnd = topSegIt->getStartPosition() - (hal_index_t)topSegIt->getLength() + 1;
        _endOffset = max((hal_index_t)0, botEnd - topEnd);
    }
    assert(_startOffset + _endOffset <= _bottomSegment->getLength());
    assert(inRange() == true);
}

BottomSegmentIteratorPtr BottomSegmentIterator::clone() const {
    assert(inRange() == true);
    BottomSegmentIteratorPtr botSegIt = getGenome()->getBottomSegmentIterator(getArrayIndex());
    if (_reversed) {
        botSegIt->toReverse();
    }
    botSegIt->slice(_startOffset, _endOffset);
    return botSegIt;
}

void BottomSegmentIterator::copy(const BottomSegmentIteratorPtr &botSegIt) {
    assert(botSegIt.get() != NULL);
    _bottomSegment->setArrayIndex(botSegIt->getGenome(), botSegIt->getArrayIndex());
    _startOffset = botSegIt->getStartOffset();
    _endOffset = botSegIt->getEndOffset();
    _reversed = botSegIt->getReversed();
}
