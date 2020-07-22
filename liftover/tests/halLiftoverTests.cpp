/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halApiTestSupport.h"
#include "halLiftoverTests.h"
#include "halBlockLiftover.h"
#include <cstdio>

using namespace std;
using namespace hal;

static void setupSharedAlignment(AlignmentPtr alignment) {
    Genome *root = alignment->addRootGenome("root");
    Genome *child1 = alignment->addLeafGenome("child1", "root", 1);
    Genome *leaf1 = alignment->addLeafGenome("leaf1", "root", 1);
    Genome *leaf2 = alignment->addLeafGenome("leaf2", "child1", 1);
    Genome *leaf3 = alignment->addLeafGenome("leaf3", "child1", 1);
    vector<Sequence::Info> seqVec(1);
    seqVec[0] = Sequence::Info("Sequence", 100, 0, 5);
    root->setDimensions(seqVec);
    seqVec[0] = Sequence::Info("Sequence", 100, 5, 7);
    child1->setDimensions(seqVec);
    seqVec[0] = Sequence::Info("Sequence", 100, 5, 0);
    leaf1->setDimensions(seqVec);
    seqVec[0] = Sequence::Info("Sequence", 70, 5, 0);
    leaf2->setDimensions(seqVec);
    seqVec[0] = Sequence::Info("Sequence", 100, 6, 0);
    leaf3->setDimensions(seqVec);
    root->setString("CAAAAGCTGCCTCGGCGTAGCCAGGTGTAAGCTGGTATTGTTCTTGTGCATCTGGGCACCATTCTCTTGTTCGTAAATAGGCGACGCTGTCTTTTGGCCG");
    child1->setString("CAAAAGCTGCCTCGGCGTAGCCAGGTGTAAGCTGGTATTGTTCTTGTGCATCTGGGCACCATTCTCTTGTTCGTAAATAGGCGACGCTGTCTTTTGGCCG");
    leaf1->setString("CAAAAGCTGCCTCGGCGTAGCCAGGTGTAAGCTGGTATTGTTCTTGTGCATCTGGGCACCATTCTCTTGTTCGTAAATAGGCGACGCTGTCTTTTGGCCG");
    leaf2->setString("ATGTGTATGCTTGGGTCAACTCTCTTTTCAGATCCGGGCGGTCGTCCGTAATTATGTGCCGAATCTCCAC");
    leaf3->setString("CAAAAGCTGCCTCGGCGTAGCCAGGTGTAAGCTGGTATTGTTCTTGTGCATCTGGGCACCATTCTCTTGTTCGTAAATAGGCGACGCTGTCTTTTGGCCG");
    // Abandon all hope, ye who enter here
    // Root
    BottomSegmentIteratorPtr botIt = root->getBottomSegmentIterator();
    botIt->setCoordinates(0, 20);
    botIt->bseg()->setChildIndex(0, 0);
    botIt->bseg()->setChildIndex(1, 0);
    botIt->bseg()->setChildReversed(0, true);
    botIt->bseg()->setChildReversed(1, true);
    botIt->bseg()->setTopParseIndex(NULL_INDEX);
    botIt->toRight();
    botIt->setCoordinates(20, 20);
    botIt->bseg()->setChildIndex(0, NULL_INDEX);
    botIt->bseg()->setChildIndex(1, 2);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, true);
    botIt->bseg()->setTopParseIndex(NULL_INDEX);
    botIt->toRight();
    botIt->setCoordinates(40, 20);
    botIt->bseg()->setChildIndex(0, 2);
    botIt->bseg()->setChildIndex(1, 1);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, false);
    botIt->bseg()->setTopParseIndex(NULL_INDEX);
    botIt->toRight();
    botIt->setCoordinates(60, 20);
    botIt->bseg()->setChildIndex(0, 3);
    botIt->bseg()->setChildIndex(1, NULL_INDEX);
    botIt->bseg()->setChildReversed(0, true);
    botIt->bseg()->setChildReversed(1, false);
    botIt->bseg()->setTopParseIndex(NULL_INDEX);
    botIt->toRight();
    botIt->setCoordinates(80, 20);
    botIt->bseg()->setChildIndex(0, NULL_INDEX);
    botIt->bseg()->setChildIndex(1, 4);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, false);
    botIt->bseg()->setTopParseIndex(NULL_INDEX);
    // Internal node
    TopSegmentIteratorPtr topIt = child1->getTopSegmentIterator();
    topIt->setCoordinates(0, 20);
    topIt->tseg()->setParentIndex(0);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(4);
    topIt->toRight();
    topIt->setCoordinates(20, 20);
    topIt->tseg()->setParentIndex(NULL_INDEX);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(40, 20);
    topIt->tseg()->setParentIndex(2);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(60, 20);
    topIt->tseg()->setParentIndex(3);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(80, 20);
    topIt->tseg()->setParentIndex(0);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(0);
    botIt = child1->getBottomSegmentIterator();
    botIt->setCoordinates(0, 20);
    botIt->bseg()->setChildIndex(0, 0);
    botIt->bseg()->setChildIndex(1, NULL_INDEX);
    botIt->bseg()->setChildReversed(0, true);
    botIt->bseg()->setChildReversed(1, false);
    botIt->toRight();
    botIt->setCoordinates(20, 10);
    botIt->bseg()->setChildIndex(0, NULL_INDEX);
    botIt->bseg()->setChildIndex(1, 0);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, true);
    botIt->toRight();
    botIt->setCoordinates(30, 5);
    botIt->bseg()->setChildIndex(0, 1);
    botIt->bseg()->setChildIndex(1, NULL_INDEX);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, false);
    botIt->toRight();
    botIt->setCoordinates(35, 15);
    botIt->bseg()->setChildIndex(0, NULL_INDEX);
    botIt->bseg()->setChildIndex(1, 2);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, false);
    botIt->toRight();
    botIt->setCoordinates(50, 20);
    botIt->bseg()->setChildIndex(0, 4);
    botIt->bseg()->setChildIndex(1, 1);
    botIt->bseg()->setChildReversed(0, true);
    botIt->bseg()->setChildReversed(1, true);
    botIt->toRight();
    botIt->setCoordinates(70, 20);
    botIt->bseg()->setChildIndex(0, 3);
    botIt->bseg()->setChildIndex(1, 3);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, true);
    botIt->toRight();
    botIt->setCoordinates(90, 10);
    botIt->bseg()->setChildIndex(0, NULL_INDEX);
    botIt->bseg()->setChildIndex(1, 4);
    botIt->bseg()->setChildReversed(0, false);
    botIt->bseg()->setChildReversed(1, false);
    // Leaf node 1 (sibling to internal)
    topIt = leaf1->getTopSegmentIterator();
    topIt->tseg()->setCoordinates(0, 20);
    topIt->tseg()->setParentIndex(0);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(20, 20);
    topIt->tseg()->setParentIndex(2);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(40, 20);
    topIt->tseg()->setParentIndex(1);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(60, 20);
    topIt->tseg()->setParentIndex(NULL_INDEX);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(80, 20);
    topIt->tseg()->setParentIndex(4);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);

    // Leaf node 2
    topIt = leaf2->getTopSegmentIterator();
    topIt->tseg()->setCoordinates(0, 20);
    topIt->tseg()->setParentIndex(0);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(20, 5);
    topIt->tseg()->setParentIndex(2);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(2);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(25, 5);
    topIt->tseg()->setParentIndex(2);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(1);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(30, 20);
    topIt->tseg()->setParentIndex(5);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(50, 20);
    topIt->tseg()->setParentIndex(4);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    // Leaf node 3
    topIt = leaf3->getTopSegmentIterator();
    topIt->tseg()->setCoordinates(0, 10);
    topIt->tseg()->setParentIndex(1);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(10, 20);
    topIt->tseg()->setParentIndex(4);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(30, 15);
    topIt->tseg()->setParentIndex(3);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(45, 20);
    topIt->tseg()->setParentIndex(5);
    topIt->tseg()->setParentReversed(true);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(65, 10);
    topIt->tseg()->setParentIndex(6);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    topIt->toRight();
    topIt->tseg()->setCoordinates(75, 25);
    topIt->tseg()->setParentIndex(NULL_INDEX);
    topIt->tseg()->setParentReversed(false);
    topIt->tseg()->setNextParalogyIndex(NULL_INDEX);
    topIt->tseg()->setBottomParseIndex(NULL_INDEX);
    root->fixParseInfo();
    child1->fixParseInfo();
    leaf1->fixParseInfo();
    leaf2->fixParseInfo();
    leaf3->fixParseInfo();
    alignment->closeGenome(root);
    alignment->closeGenome(child1);
    alignment->closeGenome(leaf1);
    alignment->closeGenome(leaf2);
    alignment->closeGenome(leaf3);
}

void BedLiftoverTest::liftAndCheck(AlignmentConstPtr alignment,
                                   const Genome *srcGenome,
                                   const Genome *tgtGenome,
                                   const string& inBed,
                                   const string& expectBed,
                                   bool outPSL, bool outPSLWithName) {
    BlockLiftover liftover;
    stringstream bedFile(inBed);
    stringstream outStream;
    liftover.convert(alignment, srcGenome, &bedFile, tgtGenome, &outStream,
                     0, true, outPSL, outPSLWithName);
    if (outStream.str() != expectBed) {
        cerr << "Got: " << endl << outStream.str() << endl;
        cerr << "Expected: " << endl << expectBed << endl;
    }
    CuAssertTrue(_testCase, outStream.str() == expectBed);
}

void BedLiftoverTest::testOneBranchLifts(AlignmentConstPtr alignment) {
    BlockLiftover liftover;
    const Genome *root = alignment->openGenome("root");
    const Genome *child1 = alignment->openGenome("child1");
    const Genome *leaf1 = alignment->openGenome("leaf1");
    const Genome *leaf2 = alignment->openGenome("leaf2");
    const Genome *leaf3 = alignment->openGenome("leaf3");

    // whole blocks (reversed, unreversed, paralogies)
    liftAndCheck(alignment, child1, root,
                 "Sequence\t0\t20\tPARALOGY1REV\t0\t+\n"
                 "Sequence\t60\t80\tREV\t0\t+\n"
                 "Sequence\t20\t40\tINSERTION\t0\t+\n"
                 "Sequence\t80\t100\tPARALOGY2\t0\t+\n",
                 // ->
                 "Sequence\t0\t20\tPARALOGY1REV\t0\t-\n"
                 "Sequence\t60\t80\tREV\t0\t-\n"
                 "Sequence\t0\t20\tPARALOGY2\t0\t+\n");

    // segment fragments (including overlapping)
    liftAndCheck(alignment, leaf1, root,
                 "Sequence\t0\t5\tNORMALREV\t0\t+\n"
                 "Sequence\t10\t30\tOVERLAP\t0\t+\n"          // 1st half reversed
                 "Sequence\t50\t70\tOVERLAPINSERTION\t0\t+\n" // 1st half reversed
                 "Sequence\t70\t100\tOVERLAPINSERTION2\t0\t+\n",
                 //->
                 "Sequence\t15\t20\tNORMALREV\t0\t-\n"
                 "Sequence\t0\t10\tOVERLAP\t0\t-\n"
                 "Sequence\t40\t50\tOVERLAP\t0\t+\n"
                 "Sequence\t20\t30\tOVERLAPINSERTION\t0\t-\n"
                 "Sequence\t80\t100\tOVERLAPINSERTION2\t0\t+\n");

    // test lifts down a branch (root->child1)
    liftAndCheck(alignment, root, child1,
                 "Sequence\t0\t10\tPARALOGY\t0\t+\n"            // reversed in one of the duplicates
                 "Sequence\t30\t50\tOVERLAPINSERTION\t0\t+\n",
                 //->
                 "Sequence\t10\t20\tPARALOGY\t0\t-\n"
                 "Sequence\t80\t90\tPARALOGY\t0\t+\n"
                 "Sequence\t40\t50\tOVERLAPINSERTION\t0\t+\n");
    alignment->closeGenome(root);
    alignment->closeGenome(child1);
    alignment->closeGenome(leaf1);
    alignment->closeGenome(leaf2);
    alignment->closeGenome(leaf3);
}

void BedLiftoverTest::testMultiBranchLifts(AlignmentConstPtr alignment) {
    BlockLiftover liftover;
    const Genome *root = alignment->openGenome("root");
    const Genome *child1 = alignment->openGenome("child1");
    const Genome *leaf1 = alignment->openGenome("leaf1");
    const Genome *leaf2 = alignment->openGenome("leaf2");
    const Genome *leaf3 = alignment->openGenome("leaf3");

    // leaf2 -> leaf3 (up 1 then down 1)
    liftAndCheck(alignment, leaf2, leaf3,
                 "Sequence\t30\t35\tREV\t0\t+\n"
                 "Sequence\t40\t60\tOVERLAP\t0\t+\n", // 1st half reversed, 2nd half double reversed (not realistic)
                 //->
                 "Sequence\t60\t65\tREV\t0\t-\n"
                 "Sequence\t45\t55\tOVERLAP\t0\t-\n"
                 "Sequence\t10\t20\tOVERLAP\t0\t+\n");

    // root->leaf2 (down 2)
    liftAndCheck(alignment, root, leaf2,
                 "Sequence\t0\t20\tBLOCK_A\t0\t+\n"
                 "Sequence\t30\t50\tBLOCK_B\t0\t+\n",
                 //->
                 "Sequence\t0\t20\tBLOCK_A\t0\t+\n"
                 "Sequence\t40\t50\tBLOCK_A\t0\t+\n");

    // leaf3->leaf1 (up 2, down 1) with blocks
    const string case3("Sequence\t0\t10\tSEGMENT_0\t0\t+\t0\t10\t128,0,0\t1\t10\t0,\n"
                       "Sequence\t10\t30\tSEGMENT_1\t0\t+\t10\t30\t128,0,0\t1\t20\t0,\n"
                       "Sequence\t30\t45\tSEGMENT_2\t0\t+\t30\t45\t128,0,0\t1\t15\t0,\n"
                       "Sequence\t45\t65\tSEGMENT_3\t0\t+\t45\t65\t128,0,0\t1\t20\t0,\n"
                       "Sequence\t65\t75\tSEGMENT_4\t0\t+\t65\t75\t128,0,0\t1\t10\t0,\n"
                       "Sequence\t75\t100\tSEGMENT_5\t0\t+\t75\t100\t128,0,0\t1\t25\t0,\n");
    liftAndCheck(alignment, leaf3, leaf1,
                 case3,
                 //->
                 "Sequence\t30\t40\tSEGMENT_1\t0\t-\t30\t40\t128,0,0\t1\t10\t0\n"
                 "Sequence\t20\t30\tSEGMENT_2\t0\t+\t20\t30\t128,0,0\t1\t10\t0\n"
                 "Sequence\t10\t20\tSEGMENT_3\t0\t+\t10\t20\t128,0,0\t1\t10\t0\n"
                 "Sequence\t0\t10\tSEGMENT_4\t0\t-\t0\t10\t128,0,0\t1\t10\t0\n");
    liftAndCheck(alignment, leaf3, leaf1,
                 case3,
                 //->
                 "2\t8\t0\t0\t0\t0\t0\t0\t+-\tSequence\t100\t20\t30\tSequence\t100\t30\t40\t1\t10,\t20,\t60,\n"
                 "2\t8\t0\t0\t0\t0\t0\t0\t++\tSequence\t100\t35\t45\tSequence\t100\t20\t30\t1\t10,\t35,\t20,\n"
                 "3\t7\t0\t0\t0\t0\t0\t0\t++\tSequence\t100\t45\t55\tSequence\t100\t10\t20\t1\t10,\t45,\t10,\n"
                 "3\t7\t0\t0\t0\t0\t0\t0\t+-\tSequence\t100\t65\t75\tSequence\t100\t0\t10\t1\t10,\t65,\t90,\n",
                 true);
    liftAndCheck(alignment, leaf3, leaf1,
                 case3,
                 //->
                 "SEGMENT_1\t2\t8\t0\t0\t0\t0\t0\t0\t+-\tSequence\t100\t20\t30\tSequence\t100\t30\t40\t1\t10,\t20,\t60,\n"
                 "SEGMENT_2\t2\t8\t0\t0\t0\t0\t0\t0\t++\tSequence\t100\t35\t45\tSequence\t100\t20\t30\t1\t10,\t35,\t20,\n"
                 "SEGMENT_3\t3\t7\t0\t0\t0\t0\t0\t0\t++\tSequence\t100\t45\t55\tSequence\t100\t10\t20\t1\t10,\t45,\t10,\n"
                 "SEGMENT_4\t3\t7\t0\t0\t0\t0\t0\t0\t+-\tSequence\t100\t65\t75\tSequence\t100\t0\t10\t1\t10,\t65,\t90,\n",
                 true, true);
    alignment->closeGenome(root);
    alignment->closeGenome(child1);
    alignment->closeGenome(leaf1);
    alignment->closeGenome(leaf2);
    alignment->closeGenome(leaf3);
}

void BedLiftoverTest::createCallBack(AlignmentPtr alignment) {
    setupSharedAlignment(alignment);
}

void BedLiftoverTest::checkCallBack(AlignmentConstPtr alignment) {
    testOneBranchLifts(alignment);
    testMultiBranchLifts(alignment);
}

/* FIXME: what is this?
// Makes assumptions about wig output:
// 1. input wig step = output wig step
// 2. only fixedStep is used
void WiggleLiftoverTest::checkWigHasEntry(vector<string> &wig,
                                          vector<string> &entry)
{
  string header = entry[0];

}
*/

void WiggleLiftoverTest::testOneBranchLifts(AlignmentConstPtr alignment) {
#if 0 // FIXME: implement
    const Genome *root = alignment->openGenome("root");
    const Genome *child1 = alignment->openGenome("child1");
    const Genome *leaf1 = alignment->openGenome("leaf1");
    const Genome *leaf2 = alignment->openGenome("leaf2");
    const Genome *leaf3 = alignment->openGenome("leaf3");

    // whole blocks (reversed, unreversed, paralogies)
    const string case1("fixedStep chrom=Sequence start=0 step=5\n"
                       "1\n2\n3\n4\n"
                       "fixedStep chrom=Sequence start=60 step=5\n"
                       "5\n6\n7\n8\n"
                       "fixedStep chrom=Sequence start=20 step=5\n"
                       "9\n10\n11\n12\n"
                       "fixedStep chrom=Sequence start=80 step=5\n"
                       "13\n14\n15\n16\n");
    stringstream wigFile(case1);
    stringstream outStream;
    liftover.convert(alignment, child1, &wigFile, root, &outStream);
    CuAssertTrue(_testCase, outStream.str() ==
                 "Sequence\t0\t20\tPARALOGY1REV\t0\t-\n"
                 "Sequence\t60\t80\tREV\t0\t-\n"
                 "Sequence\t0\t20\tPARALOGY2\t0\t+\n");
    alignment->closeGenome(root);
    alignment->closeGenome(child1);
    alignment->closeGenome(leaf1);
    alignment->closeGenome(leaf2);
    alignment->closeGenome(leaf3);
#endif
}

void WiggleLiftoverTest::testMultiBranchLifts(AlignmentConstPtr alignment) {
}

void WiggleLiftoverTest::createCallBack(AlignmentPtr alignment) {
    setupSharedAlignment(alignment);
}

void WiggleLiftoverTest::checkCallBack(AlignmentConstPtr alignment) {
    testOneBranchLifts(alignment);
    testMultiBranchLifts(alignment);
}

void halBedLiftoverTest(CuTest *testCase) {
    try {
        BedLiftoverTest tester;
        tester.check(testCase);
    } catch (...) {
        CuAssertTrue(testCase, false);
    }
}

void halWiggleLiftoverTest(CuTest *testCase) {
    try {
        WiggleLiftoverTest tester;
        tester.check(testCase);
    } catch (...) {
        CuAssertTrue(testCase, false);
    }
}

CuSuite *halLiftoverTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halBedLiftoverTest);
    SUITE_ADD_TEST(suite, halWiggleLiftoverTest);
    return suite;
}

int halLiftoverRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite *suite = CuSuiteNew();
    CuSuiteAddSuite(suite, halLiftoverTestSuite());
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
    return halLiftoverRunAllTests();
}
