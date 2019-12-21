/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halApiTestSupport.h"
#include "halSegmentTestSupport.h"
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

struct RearrangementInsertionTest : public AlignmentTest {
    void createCallBack(Alignment *alignment) {
        size_t numSequences = 3;
        size_t numSegmentsPerSequence = 10;
        size_t segmentLength = 50;

        addIdenticalParentChild(alignment, numSequences, numSegmentsPerSequence, segmentLength);

        Genome *parent = alignment->openGenome("parent");
        Genome *child = alignment->openGenome("child");

        BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();

        // insertion smaller than gap threshold
        bi->toRight();
        makeInsertion(bi);

        // stagger insertions to prevent gapped iterators from being larger
        // than a segment
        size_t count = 0;
        for (bi = parent->getBottomSegmentIterator(); not bi->atEnd(); bi->toRight()) {
            if (bi->bseg()->hasChild(0)) {
                bi->getBottomSegment()->setChildReversed(0, count % 2);
                TopSegmentIteratorPtr ti = child->getTopSegmentIterator();
                ti->toChild(bi, 0);
                ti->getTopSegment()->setParentReversed(count % 2);
                ++count;
            }
        }

        // insertion larger than gap threshold but that contains gaps
    }

    void checkCallBack(const Alignment *alignment) {
        BottomSegmentIteratorPtr bi;
        TopSegmentIteratorPtr ti;

        const Genome *child = alignment->openGenome("child");

        RearrangementPtr r = child->getRearrangement(0, 10, 1);
        do {
            hal_index_t leftIdx = r->getLeftBreakpoint()->getTopSegment()->getArrayIndex();
            if (leftIdx == 1) {
                CuAssertTrue(_testCase, r->getID() == Rearrangement::Insertion);
            } else if (leftIdx == 2) {
                // side effect of makeInsertion causes a deletion right next door
                CuAssertTrue(_testCase, r->getID() == Rearrangement::Deletion);
            } else {
                CuAssertTrue(_testCase, r->getID() != Rearrangement::Insertion && r->getID() != Rearrangement::Deletion);
            }
        } while (r->identifyNext() == true);
    }
};

struct RearrangementSimpleInversionTest : public AlignmentTest {
    void createCallBack(Alignment *alignment) {
        size_t numSequences = 3;
        size_t numSegmentsPerSequence = 10;
        size_t segmentLength = 50;

        addIdenticalParentChild(alignment, numSequences, numSegmentsPerSequence, segmentLength);

        Genome *parent = alignment->openGenome("parent");
        Genome *child = alignment->openGenome("child");

        // inversion on 1st segment (interior)
        BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator(1);
        TopSegmentIteratorPtr ti = child->getTopSegmentIterator(1);
        bi->getBottomSegment()->setChildReversed(0, true);
        ti->getTopSegment()->setParentReversed(true);

        // inversion on last segment of first sequence
        bi = parent->getBottomSegmentIterator(numSegmentsPerSequence - 1);
        ti = child->getTopSegmentIterator(numSegmentsPerSequence - 1);
        bi->getBottomSegment()->setChildReversed(0, true);
        ti->getTopSegment()->setParentReversed(true);

        // inversion on first segment of 3rd sequence
        bi = parent->getBottomSegmentIterator(numSegmentsPerSequence * 2);
        ti = child->getTopSegmentIterator(numSegmentsPerSequence * 2);
        bi->getBottomSegment()->setChildReversed(0, true);
        ti->getTopSegment()->setParentReversed(true);
    }

    void checkCallBack(const Alignment *alignment) {
        BottomSegmentIteratorPtr bi;
        TopSegmentIteratorPtr ti;

        const Genome *child = alignment->openGenome("child");

        size_t numSegmentsPerSequence = child->getSequenceIterator()->getSequence()->getNumTopSegments();

        RearrangementPtr r = child->getRearrangement(0, 10, 1);

        r = child->getRearrangement(0, 10, 1);
        do {
            hal_index_t leftIdx = r->getLeftBreakpoint()->getTopSegment()->getArrayIndex();
            if (leftIdx == 1 || leftIdx == (hal_index_t)(numSegmentsPerSequence - 1) ||
                leftIdx == (hal_index_t)(numSegmentsPerSequence * 2)) {
                CuAssertTrue(_testCase, r->getID() == Rearrangement::Inversion);
            } else {
                CuAssertTrue(_testCase, r->getID() != Rearrangement::Inversion);
            }
        } while (r->identifyNext() == true);
    }
};

struct RearrangementGappedInversionTest : public AlignmentTest {
    void createCallBack(Alignment *alignment) {
        size_t numSequences = 3;
        size_t numSegmentsPerSequence = 10;
        size_t segmentLength = 5;

        addIdenticalParentChild(alignment, numSequences, numSegmentsPerSequence, segmentLength);

        Genome *parent = alignment->openGenome("parent");
        Genome *child = alignment->openGenome("child");

        // 4-segment inversion.
        // parent has gap-deletions at 2nd and 5th posstions
        // child has gap-insertions positions 3 and 5
        BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator(1);
        TopSegmentIteratorPtr ti = child->getTopSegmentIterator(1);
        bi->getBottomSegment()->setChildReversed(0, true);
        bi->getBottomSegment()->setChildIndex(0, 6);
        ti->getTopSegment()->setParentReversed(true);
        ti->getTopSegment()->setParentIndex(6);

        bi = parent->getBottomSegmentIterator(2);
        ti = child->getTopSegmentIterator(2);
        bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
        ti->getTopSegment()->setParentReversed(true);
        ti->getTopSegment()->setParentIndex(4);

        bi = parent->getBottomSegmentIterator(3);
        ti = child->getTopSegmentIterator(3);
        bi->getBottomSegment()->setChildReversed(0, true);
        bi->getBottomSegment()->setChildIndex(0, 4);
        ti->getTopSegment()->setParentIndex(NULL_INDEX);

        bi = parent->getBottomSegmentIterator(4);
        ti = child->getTopSegmentIterator(4);
        bi->getBottomSegment()->setChildReversed(0, true);
        bi->getBottomSegment()->setChildIndex(0, 2);
        ti->getTopSegment()->setParentIndex(3);
        ti->getTopSegment()->setParentReversed(true);

        bi = parent->getBottomSegmentIterator(5);
        ti = child->getTopSegmentIterator(5);
        bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
        ti->getTopSegment()->setParentIndex(NULL_INDEX);

        bi = parent->getBottomSegmentIterator(6);
        ti = child->getTopSegmentIterator(6);
        bi->getBottomSegment()->setChildReversed(0, true);
        bi->getBottomSegment()->setChildIndex(0, 1);
        ti->getTopSegment()->setParentIndex(1);
        ti->getTopSegment()->setParentReversed(true);
    }

    void checkCallBack(const Alignment *alignment) {
        BottomSegmentIteratorPtr bi;
        TopSegmentIteratorPtr ti;

        const Genome *child = alignment->openGenome("child");

        RearrangementPtr r = child->getRearrangement(0, 10, 1);
        do {
            hal_index_t leftIdx = r->getLeftBreakpoint()->getTopSegment()->getArrayIndex();
            if (leftIdx == 1) {
                CuAssertTrue(_testCase, r->getID() == Rearrangement::Inversion);
            } else {
                CuAssertTrue(_testCase, r->getID() != Rearrangement::Inversion);
            }
        } while (r->identifyNext() == true);
    }
};

static void halRearrangementInsertionTest(CuTest *testCase) {
    RearrangementInsertionTest tester;
    tester.check(testCase);
}

static void halRearrangementSimpleInversionTest(CuTest *testCase) {
    RearrangementSimpleInversionTest tester;
    tester.check(testCase);
}

static void halRearrangementGappedInversionTest(CuTest *testCase) {
    RearrangementGappedInversionTest tester;
    tester.check(testCase);
}

static CuSuite *halRearrangementTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halRearrangementInsertionTest);
    SUITE_ADD_TEST(suite, halRearrangementSimpleInversionTest);
    SUITE_ADD_TEST(suite, halRearrangementGappedInversionTest);
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halRearrangementTestSuite());
}
