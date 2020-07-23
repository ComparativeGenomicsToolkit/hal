/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halApiTestSupport.h"
#include <iostream>
#include <string>
#include "halSequence.h"
#include "halGenome.h"
#include "halSequenceIterator.h"
#include "halTopSegmentIterator.h"
#include "halBottomSegmentIterator.h"

extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

struct SequenceCreateTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        hal_size_t alignmentSize = alignment->getNumGenomes();
        CuAssertTrue(_testCase, alignmentSize == 0);

        Genome *ancGenome = alignment->addRootGenome("AncGenome", 0);

        size_t numSequences = 1000;
        vector<Sequence::Info> seqVec;
        for (size_t i = 0; i < numSequences; ++i) {
            hal_size_t len = 1 + i * 5 + i;
            seqVec.push_back(Sequence::Info("sequence" + std::to_string(i), len, i, i * 2));
        }
        ancGenome->setDimensions(seqVec);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        const Genome *ancGenome = alignment->openGenome("AncGenome");

        hal_size_t numSequences = ancGenome->getNumSequences();
        CuAssertTrue(_testCase, numSequences = 1000);

        hal_size_t numTopSegments = 0;
        hal_size_t numBottomSegments = 0;
        hal_size_t totalLength = 0;
        hal_size_t lastStart = 0;
        hal_size_t lastLength = 0;

        for (SequenceIteratorPtr seqIt = ancGenome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            hal_size_t i = (hal_size_t)seqIt->getSequence()->getArrayIndex();
            hal_size_t len = 1 + i * 5 + i;
            string name = "sequence" + std::to_string(i);
            const Sequence *seq = seqIt->getSequence();
            CuAssertTrue(_testCase, seq->getName() == name);
            CuAssertTrue(_testCase, seq->getSequenceLength() == len);
            CuAssertTrue(_testCase, seq->getNumTopSegments() == i);
            CuAssertTrue(_testCase, seq->getNumBottomSegments() == i * 2);
            const Genome *gen = seq->getGenome();
            CuAssertTrue(_testCase, gen->getName() == "AncGenome");

            numTopSegments += seq->getNumTopSegments();
            numBottomSegments += seq->getNumBottomSegments();
            totalLength += seq->getSequenceLength();

            if (i == 0) {
                CuAssertTrue(_testCase, seq->getStartPosition() == 0);
            } else {
                CuAssertTrue(_testCase, seq->getStartPosition() - lastStart == lastLength);
            }

            lastStart = seq->getStartPosition();
            lastLength = seq->getSequenceLength();
        }

        const Sequence *seq = ancGenome->getSequence("sequence555");
        CuAssertTrue(_testCase, seq->getName() == "sequence555");
        seq = ancGenome->getSequenceBySite(0);
        CuAssertTrue(_testCase, seq->getName() == "sequence0");
        seq = ancGenome->getSequenceBySite(45);
        CuAssertTrue(_testCase, seq->getName() == "sequence4");

        CuAssertTrue(_testCase, ancGenome->getSequenceLength() == totalLength);
        CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == numTopSegments);
        CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == numBottomSegments);
    }
};

struct SequenceIteratorTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        hal_size_t alignmentSize = alignment->getNumGenomes();
        CuAssertTrue(_testCase, alignmentSize == 0);

        Genome *ancGenome = alignment->addRootGenome("AncGenome", 0);

        size_t numSequences = 1000;
        vector<Sequence::Info> seqVec;
        for (size_t i = 0; i < numSequences; ++i) {
            size_t len = 100;
            string name = "sequence" + std::to_string(i);
            seqVec.push_back(Sequence::Info(name, i * len, i ? len : 0, i ? len : 0));
        }
        ancGenome->setDimensions(seqVec);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        const Genome *ancGenome = alignment->openGenome("AncGenome");

        hal_size_t numSequences = ancGenome->getNumSequences();
        CuAssertTrue(_testCase, numSequences = 1000);

        for (SequenceIteratorPtr seqIt = ancGenome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *seq = seqIt->getSequence();
            hal_size_t i = seq->getArrayIndex();

            TopSegmentIteratorPtr tsIt = seq->getTopSegmentIterator();
            hal_size_t numTopSegments = seq->getNumTopSegments();
            for (hal_size_t j = 0; j < numTopSegments; ++j) {
                TopSegmentIteratorPtr gtsIt = ancGenome->getTopSegmentIterator((i - 1) * 100 + j);
                const TopSegment *gsTopSegment = gtsIt->getTopSegment();
                const TopSegment *sqTopSegment = tsIt->getTopSegment();

                CuAssertTrue(_testCase, gsTopSegment->getArrayIndex() == sqTopSegment->getArrayIndex());
                tsIt->toRight();
            }

            BottomSegmentIteratorPtr bsIt = seq->getBottomSegmentIterator();
            hal_size_t numBottomSegments = seq->getNumBottomSegments();
            for (hal_size_t j = 0; j < numBottomSegments; ++j) {
                BottomSegmentIteratorPtr gbsIt = ancGenome->getBottomSegmentIterator((i - 1) * 100 + j);
                const BottomSegment *gsBottomSegment = gbsIt->getBottomSegment();
                const BottomSegment *sqBottomSegment = bsIt->getBottomSegment();

                CuAssertTrue(_testCase, gsBottomSegment->getArrayIndex() == sqBottomSegment->getArrayIndex());
                bsIt->toRight();
            }
        }
    }
};

struct SequenceUpdateTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        hal_size_t alignmentSize = alignment->getNumGenomes();
        CuAssertTrue(_testCase, alignmentSize == 0);

        Genome *ancGenome = alignment->addRootGenome("AncGenome", 0);

        size_t numSequences = 1000;
        vector<Sequence::Info> seqVec;
        for (size_t i = 0; i < numSequences; ++i) {
            hal_size_t len = 1 + i * 5 + i;
            string name = "sequence" + std::to_string(i);
            seqVec.push_back(Sequence::Info(name, len, i, i * 2));
        }
        ancGenome->setDimensions(seqVec);
        alignment->closeGenome(ancGenome);
        ancGenome = alignment->openGenome("AncGenome");

        vector<Sequence::UpdateInfo> updateVec;
        for (size_t i = 0; i < numSequences / 2; ++i) {
            const Sequence *sequence = ancGenome->getSequence("sequence" + std::to_string(i));
            updateVec.push_back(Sequence::UpdateInfo(sequence->getName(), i * 7));
        }
        ancGenome->updateTopDimensions(updateVec);

        updateVec.clear();
        for (size_t i = 0; i < numSequences / 3; ++i) {
            const Sequence *sequence = ancGenome->getSequence("sequence" + std::to_string(i));
            updateVec.push_back(Sequence::UpdateInfo(sequence->getName(), i * 5));
        }
        ancGenome->updateBottomDimensions(updateVec);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        const Genome *ancGenome = alignment->openGenome("AncGenome");

        hal_size_t numSequences = ancGenome->getNumSequences();
        CuAssertTrue(_testCase, numSequences = 1000);

        hal_size_t numTopSegments = 0;
        hal_size_t numBottomSegments = 0;
        hal_size_t totalLength = 0;
        hal_size_t lastStart = 0;
        hal_size_t lastLength = 0;

        for (SequenceIteratorPtr seqIt = ancGenome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            hal_size_t i = (hal_size_t)sequence->getArrayIndex();
            hal_size_t len = 1 + i * 5 + i;
            hal_size_t numTop = i < numSequences / 2 ? i * 7 : i;
            hal_size_t numBottom = i < numSequences / 3 ? i * 5 : i * 2;
            string name = "sequence" + std::to_string(i);
            const Sequence *seq = seqIt->getSequence();
            CuAssertTrue(_testCase, seq->getName() == name);
            CuAssertTrue(_testCase, seq->getSequenceLength() == len);
            CuAssertTrue(_testCase, seq->getNumTopSegments() == numTop);
            CuAssertTrue(_testCase, seq->getNumBottomSegments() == numBottom);
            const Genome *gen = seq->getGenome();
            CuAssertTrue(_testCase, gen->getName() == "AncGenome");

            numTopSegments += seq->getNumTopSegments();
            numBottomSegments += seq->getNumBottomSegments();
            totalLength += seq->getSequenceLength();

            if (i == 0) {
                CuAssertTrue(_testCase, seq->getStartPosition() == 0);
            } else {
                CuAssertTrue(_testCase, seq->getStartPosition() - lastStart == lastLength);
            }

            lastStart = seq->getStartPosition();
            lastLength = seq->getSequenceLength();
        }

        const Sequence *seq = ancGenome->getSequence("sequence555");
        CuAssertTrue(_testCase, seq->getName() == "sequence555");
        seq = ancGenome->getSequenceBySite(0);
        CuAssertTrue(_testCase, seq->getName() == "sequence0");
        seq = ancGenome->getSequenceBySite(45);
        CuAssertTrue(_testCase, seq->getName() == "sequence4");

        CuAssertTrue(_testCase, ancGenome->getSequenceLength() == totalLength);
        CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == numTopSegments);
        CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == numBottomSegments);
    }
};

struct SequenceRenameTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        Genome *ancGenome = alignment->addRootGenome("AncGenome", 0);

        vector<Sequence::Info> seqVec;
        seqVec.push_back(Sequence::Info("seq1", 40, 20, 20));
        seqVec.push_back(Sequence::Info("seq2", 50, 20, 20));
        seqVec.push_back(Sequence::Info("anotherSeq", 60, 20, 20));
        ancGenome->setDimensions(seqVec);
        alignment->closeGenome(ancGenome);

        ancGenome = alignment->openGenome("AncGenome");
        ancGenome->getSequence("seq1")->setName("foobar");
        ancGenome->getSequence("seq2")->setName("a_really_really_super_long_name");
        ancGenome->getSequence("anotherSeq")->setName("short");
        ancGenome->getSequence("short")->setName("again");
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        const Genome *ancGenome = alignment->openGenome("AncGenome");

        CuAssertTrue(_testCase, ancGenome->getNumSequences() == 3);
        CuAssertTrue(_testCase, ancGenome->getSequence("seq1") == NULL);
        CuAssertTrue(_testCase, ancGenome->getSequence("seq2") == NULL);
        CuAssertTrue(_testCase, ancGenome->getSequence("anotherSeq") == NULL);
        CuAssertTrue(_testCase, ancGenome->getSequence("short") == NULL);

        CuAssertTrue(_testCase, ancGenome->getSequence("foobar") != NULL);
        CuAssertTrue(_testCase, ancGenome->getSequence("a_really_really_super_long_name") != NULL);
        CuAssertTrue(_testCase, ancGenome->getSequence("again") != NULL);

        CuAssertTrue(_testCase, ancGenome->getSequence("foobar")->getSequenceLength() == 40);
        CuAssertTrue(_testCase, ancGenome->getSequence("a_really_really_super_long_name")->getSequenceLength() == 50);
        CuAssertTrue(_testCase, ancGenome->getSequence("again")->getSequenceLength() == 60);
    }
};

static void halSequenceCreateTest(CuTest *testCase) {
    SequenceCreateTest tester;
    tester.check(testCase);
}

static void halSequenceIteratorTest(CuTest *testCase) {
    SequenceIteratorTest tester;
    tester.check(testCase);
}

static void halSequenceUpdateTest(CuTest *testCase) {
    SequenceUpdateTest tester;
    tester.check(testCase);
}

static void halSequenceRenameTest(CuTest *testCase) {
    SequenceRenameTest tester;
    tester.check(testCase);
}

static CuSuite *halSequenceTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halSequenceCreateTest);
    SUITE_ADD_TEST(suite, halSequenceIteratorTest);
    SUITE_ADD_TEST(suite, halSequenceUpdateTest);
    SUITE_ADD_TEST(suite, halSequenceRenameTest);
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halSequenceTestSuite());
}
