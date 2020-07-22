/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEGMENTTESTSUPPORT_H
#define _HALSEGMENTTESTSUPPORT_H

#include "hal.h"
#include <vector>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace hal;

/* just create a bunch of garbage data.  we don't care about logical
 * consistency for this test, just whether or not it's read and written
 * properly.
 */
struct TopSegmentStruct {
    hal_size_t _length;
    hal_index_t _startPosition;
    hal_index_t _nextParalogyIndex;
    hal_index_t _parentIndex;
    bool _parentReversed;
    hal_index_t _arrayIndex;
    hal_index_t _bottomParseIndex;

    void setRandom() {
        _length = rand();
        _startPosition = rand();
        _nextParalogyIndex = rand();
        _parentIndex = rand();
        _arrayIndex = rand();
        _bottomParseIndex = rand();
    }

    void set(hal_index_t startPosition, hal_size_t length, hal_index_t parentIndex = NULL_INDEX, bool parentReversed = false,
             hal_index_t bottomParseIndex = NULL_INDEX, hal_index_t nextParalogyIndex = NULL_INDEX) {
        _startPosition = startPosition;
        _length = length;
        _parentIndex = parentIndex;
        _parentReversed = parentReversed;
        _bottomParseIndex = bottomParseIndex;
        _nextParalogyIndex = nextParalogyIndex;
    }

    void applyTo(TopSegmentIteratorPtr it) const {
        TopSegment *seg = it->getTopSegment();
        seg->setCoordinates(_startPosition, _length);
        seg->setNextParalogyIndex(_nextParalogyIndex);
        seg->setParentIndex(_parentIndex);
        seg->setParentReversed(_parentReversed);
        seg->setBottomParseIndex(_bottomParseIndex);
    }
    
    void compareTo(TopSegmentIteratorPtr it, CuTest *testCase) const {
        const TopSegment *seg = it->getTopSegment();
        CuAssertTrue(testCase, _length == seg->getLength());
        CuAssertTrue(testCase, _startPosition == seg->getStartPosition());
        CuAssertTrue(testCase, _nextParalogyIndex == seg->getNextParalogyIndex());
        CuAssertTrue(testCase, _parentIndex == seg->getParentIndex());
        CuAssertTrue(testCase, _bottomParseIndex == seg->getBottomParseIndex());
    }
};

struct BottomSegmentStruct {
    hal_size_t _length;
    hal_index_t _startPosition;
    std::vector<std::pair<hal_index_t, bool>> _children;
    hal_index_t _arrayIndex;
    hal_index_t _topParseIndex;

    // just create a bunch of garbage data.  we don't care
    // about logical consistency for this test, just whether or not
    // it's read and written properly.
    void setRandom(hal_size_t numChildren) {
        _length = rand();
        _startPosition = rand();
        _arrayIndex = rand();
        _topParseIndex = rand();
        _children.clear();
        for (hal_size_t i = 0; i < numChildren; ++i) {
            pair<hal_index_t, bool> child;
            child.first = (bool)rand() % 2;
            child.second = rand();
            _children.push_back(child);
        }
    }

    void set(hal_index_t startPosition, hal_size_t length, hal_index_t topParseIndex = NULL_INDEX) {
        _startPosition = startPosition;
        _length = length;
        _topParseIndex = topParseIndex;
    }

    void applyTo(BottomSegmentIteratorPtr it) const {
        BottomSegment *seg = it->getBottomSegment();
        seg->setCoordinates(_startPosition, _length);
        seg->setTopParseIndex(_topParseIndex);
        for (hal_size_t i = 0; i < _children.size(); ++i) {
            seg->setChildIndex(i, _children[i].first);
            seg->setChildReversed(i, _children[i].second);
        }
    }

    void compareTo(BottomSegmentIteratorPtr it, CuTest *testCase) const {
        const BottomSegment *seg = it->getBottomSegment();
        CuAssertTrue(testCase, _length == seg->getLength());
        CuAssertTrue(testCase, _startPosition == seg->getStartPosition());
        CuAssertTrue(testCase, _topParseIndex == seg->getTopParseIndex());
        CuAssertTrue(testCase, _children.size() == seg->getNumChildren());
        for (hal_size_t i = 0; i < _children.size(); ++i) {
            CuAssertTrue(testCase, _children[i].first == seg->getChildIndex(i));
            CuAssertTrue(testCase, _children[i].second == seg->getChildReversed(i));
        }
    }
};

static inline void addIdenticalParentChild(AlignmentPtr alignment, size_t numSequences, size_t numSegmentsPerSequence, size_t segmentLength) {
    vector<Sequence::Info> seqVec(numSequences);

    BottomSegmentIteratorPtr bi;
    BottomSegmentStruct bs;
    TopSegmentIteratorPtr ti;
    TopSegmentStruct ts;

    Genome *parent = alignment->addRootGenome("parent");
    Genome *child = alignment->addLeafGenome("child", "parent", 1);

    for (size_t i = 0; i < numSequences; ++i) {
        seqVec[i] = Sequence::Info("Sequence" + std::to_string(i), segmentLength * numSegmentsPerSequence,
                                   numSegmentsPerSequence, numSegmentsPerSequence);
    }
    parent->setDimensions(seqVec);
    child->setDimensions(seqVec);

    for (bi = parent->getBottomSegmentIterator(); not bi->atEnd(); bi->toRight()) {
        bs.set(bi->getBottomSegment()->getArrayIndex() * segmentLength, segmentLength);
        bs._children.clear();
        bs._children.push_back(pair<hal_size_t, bool>(bi->getBottomSegment()->getArrayIndex(), false));
        bs.applyTo(bi);
    }

    for (ti = child->getTopSegmentIterator(); not ti->atEnd(); ti->toRight()) {
        ts.set(ti->getTopSegment()->getArrayIndex() * segmentLength, segmentLength, ti->getTopSegment()->getArrayIndex());
        ts.applyTo(ti);
    }
}


// doesn't currently work on sequence endpoints
static inline void makeInsertion(BottomSegmentIteratorPtr bi) {
    assert(bi->getBottomSegment()->isLast() == false);
    TopSegmentIteratorPtr ti = bi->getBottomSegment()->getGenome()->getTopSegmentIterator();
    ti->toChild(bi, 0);
    hal_index_t pi = ti->getTopSegment()->getParentIndex();
    ti->getTopSegment()->setParentIndex(NULL_INDEX);
    ti->toRight();
    ti->getTopSegment()->setParentIndex(pi);

    hal_index_t ci = bi->getBottomSegment()->getChildIndex(0);
    bi->getBottomSegment()->setChildIndex(0, ci + 1);
    bi->toRight();
    bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
}

// designed to be used only on alignments made with
// addIdenticalParentChild()
static inline void makeInsGap(TopSegmentIteratorPtr ti) {
    Genome *genome = ti->getTopSegment()->getGenome();
    Genome *parent = genome->getParent();
    BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
    assert(ti->tseg()->hasParent() == true);
    bi->toParent(ti);
    while (not bi->atEnd()) {
        hal_index_t childIndex = bi->getBottomSegment()->getChildIndex(0);
        if (childIndex == (hal_index_t)(genome->getNumTopSegments() - 1)) {
            bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
        } else {
            bi->getBottomSegment()->setChildIndex(0, childIndex + 1);
        }
        bi->getReversed() ? bi->toLeft() : bi->toRight();
    }
    TopSegmentIteratorPtr topIt = ti->clone();
    topIt->getTopSegment()->setParentIndex(NULL_INDEX);
    topIt->toRight();
    while (not topIt->atEnd()) {
        hal_index_t parentIndex = topIt->getTopSegment()->getParentIndex();
        topIt->getTopSegment()->setParentIndex(parentIndex - 1);
        topIt->getReversed() ? topIt->toLeft() : topIt->toRight();
    }
}

static inline void makeDelGap(BottomSegmentIteratorPtr botIt) {
    Genome *parent = botIt->getBottomSegment()->getGenome();
    Genome *genome = parent->getChild(0);
    BottomSegmentIteratorPtr bi = botIt->clone();
    BottomSegmentIteratorPtr bi2 = botIt->clone();
    TopSegmentIteratorPtr ti = genome->getTopSegmentIterator();
    TopSegmentIteratorPtr ti2 = genome->getTopSegmentIterator();
    assert(bi->bseg()->hasChild(0) == true);
    hal_index_t prevIndex = bi->getBottomSegment()->getChildIndex(0);
    bool prevReversed = bi->getBottomSegment()->getChildReversed(0);
    ti->toChild(bi, 0);
    bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
    bi->toRight();
    while (not bi->atEnd()) {
        hal_index_t prevIndex2 = bi->getBottomSegment()->getChildIndex(0);
        bool prevReversed2 = bi->getBottomSegment()->getChildReversed(0);
        bi2 = bi->clone();
        assert(bi2->getReversed() == false);
        bi2->toLeft();
        bi->getBottomSegment()->setChildIndex(0, prevIndex);
        bi->getBottomSegment()->setChildReversed(0, prevReversed);
        bi->getReversed() ? bi->toLeft() : bi->toRight();
        swap(prevIndex, prevIndex2);
        swap(prevReversed, prevReversed2);
    }

    while (not ti->atEnd()) {
        hal_index_t parentIndex = ti->getTopSegment()->getParentIndex();
        if (parentIndex == (hal_index_t)(parent->getNumBottomSegments() - 1)) {
            ti->getTopSegment()->setParentIndex(NULL_INDEX);
        } else {
            ti2 = ti->clone();
            if (ti->getTopSegment()->isLast() == false) {
                ti2->getReversed() ? ti2->toLeft() : ti2->toRight();
                hal_index_t nextIndex = ti2->getTopSegment()->getParentIndex();
                bool nextRev = ti2->getTopSegment()->getParentReversed();
                ti->getTopSegment()->setParentIndex(nextIndex);
                ti->getTopSegment()->setParentReversed(nextRev);
            }
        }
        ti->getReversed() ? ti->toLeft() : ti->toRight();
    }
}

static inline void makeInversion(TopSegmentIteratorPtr ti, hal_size_t len) {
    Genome *child = ti->getTopSegment()->getGenome();
    Genome *parent = child->getParent();
    hal_index_t first = ti->getTopSegment()->getArrayIndex();
    hal_index_t last = ti->getTopSegment()->getArrayIndex() + len - 1;
    for (size_t i = 0; i < len; ++i) {
        TopSegmentIteratorPtr t = child->getTopSegmentIterator(first + i);
        BottomSegmentIteratorPtr b = parent->getBottomSegmentIterator(last - i);
        t->getTopSegment()->setParentIndex(last - i);
        t->getTopSegment()->setParentReversed(true);
        b->getBottomSegment()->setChildIndex(0, first + i);
        b->getBottomSegment()->setChildReversed(0, true);
    }
}

#endif
// Local Variables:
// mode: c++
// End:
