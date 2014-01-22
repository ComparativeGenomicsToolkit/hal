/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include "hal.h"
#include "halBlockLiftover.h"
#include "halLiftoverTests.h"

using namespace std;
using namespace hal;

void setupSharedAlignment(AlignmentPtr alignment)
{
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
  botIt->setChildIndex(0, 0);
  botIt->setChildIndex(1, 0);
  botIt->setChildReversed(0, true);
  botIt->setChildReversed(1, true);
  botIt->setTopParseIndex(NULL_INDEX);
  botIt->toRight();
  botIt->setCoordinates(20, 20);
  botIt->setChildIndex(0, NULL_INDEX);
  botIt->setChildIndex(1, 2);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, true);
  botIt->setTopParseIndex(NULL_INDEX);
  botIt->toRight();
  botIt->setCoordinates(40, 20);
  botIt->setChildIndex(0, 2);
  botIt->setChildIndex(1, 1);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);
  botIt->toRight();
  botIt->setCoordinates(60, 20);
  botIt->setChildIndex(0, 3);
  botIt->setChildIndex(1, NULL_INDEX);
  botIt->setChildReversed(0, true);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);
  botIt->toRight();
  botIt->setCoordinates(80, 20);
  botIt->setChildIndex(0, NULL_INDEX);
  botIt->setChildIndex(1, 4);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);
  // Internal node
  TopSegmentIteratorPtr topIt = child1->getTopSegmentIterator();
  topIt->setCoordinates(0, 20);
  topIt->setParentIndex(0);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(4);
  topIt->toRight();
  topIt->setCoordinates(20, 20);
  topIt->setParentIndex(NULL_INDEX);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(40, 20);
  topIt->setParentIndex(2);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(60, 20);
  topIt->setParentIndex(3);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(80, 20);
  topIt->setParentIndex(0);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(0);
  botIt = child1->getBottomSegmentIterator();
  botIt->setCoordinates(0, 20);
  botIt->setChildIndex(0, 0);
  botIt->setChildIndex(1, NULL_INDEX);
  botIt->setChildReversed(0, true);
  botIt->setChildReversed(1, false);
  botIt->toRight();
  botIt->setCoordinates(20, 10);
  botIt->setChildIndex(0, NULL_INDEX);
  botIt->setChildIndex(1, 0);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, true);
  botIt->toRight();
  botIt->setCoordinates(30, 5);
  botIt->setChildIndex(0, 1);
  botIt->setChildIndex(1, NULL_INDEX);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->toRight();
  botIt->setCoordinates(35, 15);
  botIt->setChildIndex(0, NULL_INDEX);
  botIt->setChildIndex(1, 2);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->toRight();
  botIt->setCoordinates(50, 20);
  botIt->setChildIndex(0, 4);
  botIt->setChildIndex(1, 1);
  botIt->setChildReversed(0, true);
  botIt->setChildReversed(1, true);
  botIt->toRight();
  botIt->setCoordinates(70, 20);
  botIt->setChildIndex(0, 3);
  botIt->setChildIndex(1, 3);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, true);
  botIt->toRight();
  botIt->setCoordinates(90, 10);
  botIt->setChildIndex(0, NULL_INDEX);
  botIt->setChildIndex(1, 4);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  // Leaf node 1 (sibling to internal)
  topIt = leaf1->getTopSegmentIterator();
  topIt->setCoordinates(0, 20);
  topIt->setParentIndex(0);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(20, 20);
  topIt->setParentIndex(2);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(40, 20);
  topIt->setParentIndex(1);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(60, 20);
  topIt->setParentIndex(NULL_INDEX);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(80, 20);
  topIt->setParentIndex(4);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);

  // Leaf node 2
  topIt = leaf2->getTopSegmentIterator();
  topIt->setCoordinates(0, 20);
  topIt->setParentIndex(0);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(20, 5);
  topIt->setParentIndex(2);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(2);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(25, 5);
  topIt->setParentIndex(2);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(1);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(30, 20);
  topIt->setParentIndex(5);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(50, 20);
  topIt->setParentIndex(4);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  // Leaf node 3
  topIt = leaf3->getTopSegmentIterator();
  topIt->setCoordinates(0, 10);
  topIt->setParentIndex(1);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(10, 20);
  topIt->setParentIndex(4);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(30, 15);
  topIt->setParentIndex(3);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(45, 20);
  topIt->setParentIndex(5);
  topIt->setParentReversed(true);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(65, 10);
  topIt->setParentIndex(6);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(75, 25);
  topIt->setParentIndex(NULL_INDEX);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  root->fixParseInfo();
  child1->fixParseInfo();
  leaf1->fixParseInfo();
  leaf2->fixParseInfo();
  leaf3->fixParseInfo();
}

void BedLiftoverTest::testOneBranchLifts(AlignmentConstPtr alignment)
{
  BlockLiftover liftover;
  const Genome *root = alignment->openGenome("root");
  const Genome *child1 = alignment->openGenome("child1");
  const Genome *leaf1 = alignment->openGenome("leaf1");
  const Genome *leaf2 = alignment->openGenome("leaf2");
  const Genome *leaf3 = alignment->openGenome("leaf3");

  // whole blocks (reversed, unreversed, paralogies)
  stringstream bedFile("Sequence\t0\t20\tPARALOGY1REV\t0\t+\n"
                       "Sequence\t60\t80\tREV\t0\t+\n"
                       "Sequence\t20\t40\tINSERTION\t0\t+\n"
                       "Sequence\t80\t100\tPARALOGY2\t0\t+\n");
  stringstream outStream;
  liftover.convert(alignment, child1, &bedFile, root, &outStream);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t20\tPARALOGY1REV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t60\t80\tREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t20\tPARALOGY2\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 3);

  // segment fragments (including overlapping)
  bedFile.str("Sequence\t0\t5\tNORMALREV\t0\t+\n"
              "Sequence\t10\t30\tOVERLAP\t0\t+\n" // 1st half reversed
              "Sequence\t50\t70\tOVERLAPINSERTION\t0\t+\n" // 1st half reversed
              "Sequence\t70\t100\tOVERLAPINSERTION2\t0\t+\n");
  bedFile.clear();
  outStream.str("");
  outStream.clear();
  liftover.convert(alignment, leaf1, &bedFile, root, &outStream);
  streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t15\t20\tNORMALREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t10\tOVERLAP\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t40\t50\tOVERLAP\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t20\t30\tOVERLAPINSERTION\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t80\t100\tOVERLAPINSERTION2\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 5);

  // test lifts down a branch (root->child1)
  bedFile.str("Sequence\t0\t10\tPARALOGY\t0\t+\n" // reversed in one of the duplicates
              "Sequence\t30\t50\tOVERLAPINSERTION\t0\t+\n"); // 1st half should not map
  bedFile.clear();
  outStream.str("");
  outStream.clear();
  liftover.convert(alignment, root, &bedFile, child1, &outStream);
  streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t10\t20\tPARALOGY\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t80\t90\tPARALOGY\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t40\t50\tOVERLAPINSERTION\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 3);
}

void BedLiftoverTest::testMultiBranchLifts(AlignmentConstPtr alignment)
{
  BlockLiftover liftover;
  const Genome *root = alignment->openGenome("root");
  const Genome *child1 = alignment->openGenome("child1");
  const Genome *leaf1 = alignment->openGenome("leaf1");
  const Genome *leaf2 = alignment->openGenome("leaf2");
  const Genome *leaf3 = alignment->openGenome("leaf3");

  // leaf2 -> leaf3 (up 1 then down 1)
  stringstream bedFile("Sequence\t30\t35\tREV\t0\t+\n"
                       "Sequence\t40\t60\tOVERLAP\t0\t+\n" // 1st half reversed, 2nd half double reversed (not realistic)
    );
  stringstream outStream;
  liftover.convert(alignment, leaf2, &bedFile, leaf3, &outStream);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t60\t65\tREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t45\t55\tOVERLAP\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t10\t20\tOVERLAP\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 3);

  // root->leaf2 (down 2)
  bedFile.str("Sequence\t0\t20\tBLOCK_A\t0\t+\n"
              "Sequence\t30\t50\tBLOCK_B\t0\t+\n");
  bedFile.clear();
  outStream.str("");
  outStream.clear();
  liftover.convert(alignment, root, &bedFile, leaf2, &outStream);
  streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t20\tBLOCK_A\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t40\t50\tBLOCK_A\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 2);

  // leaf3->leaf1 (up 2, down 1)
  bedFile.str("Sequence\t0\t10\tSEGMENT_0\t0\t+\n"
              "Sequence\t10\t30\tSEGMENT_1\t0\t+\n"
              "Sequence\t30\t45\tSEGMENT_2\t0\t+\n"
              "Sequence\t45\t65\tSEGMENT_3\t0\t+\n"
              "Sequence\t65\t75\tSEGMENT_4\t0\t+\n"
              "Sequence\t75\t100\tSEGMENT_5\t0\t+\n");
  bedFile.clear();
  outStream.str("");
  outStream.clear();
  liftover.convert(alignment, leaf3, &bedFile, leaf1, &outStream);
  streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t30\t40\tSEGMENT_1\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t20\t30\tSEGMENT_2\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t10\t20\tSEGMENT_3\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t10\tSEGMENT_4\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 4);
  
}

void BedLiftoverTest::createCallBack(AlignmentPtr alignment)
{
  setupSharedAlignment(alignment);
}

void BedLiftoverTest::checkCallBack(AlignmentConstPtr alignment)
{
  testOneBranchLifts(alignment);
  testMultiBranchLifts(alignment);
}

/*
// Makes assumptions about wig output:
// 1. input wig step = output wig step
// 2. only fixedStep is used
void WiggleLiftoverTest::checkWigHasEntry(vector<string> &wig,
                                          vector<string> &entry)
{
  string header = entry[0];
  
}
*/

void WiggleLiftoverTest::testOneBranchLifts(AlignmentConstPtr alignment)
{
/*  WiggleLiftover liftover;
  const Genome *root = alignment->openGenome("root");
  const Genome *child1 = alignment->openGenome("child1");
  const Genome *leaf1 = alignment->openGenome("leaf1");
  const Genome *leaf2 = alignment->openGenome("leaf2");
  const Genome *leaf3 = alignment->openGenome("leaf3");

  // whole blocks (reversed, unreversed, paralogies)
  stringstream wigFile("fixedStep chrom=Sequence start=0 step=5\n"
                       "1\n2\n3\n4\n"
                       "fixedStep chrom=Sequence start=60 step=5\n"
                       "5\n6\n7\n8\n"
                       "fixedStep chrom=Sequence start=20 step=5\n"
                       "9\n10\n11\n12\n"
                       "fixedStep chrom=Sequence start=80 step=5\n"
                       "13\n14\n15\n16\n");
  stringstream outStream;
  liftover.convert(alignment, child1, &wigFile, root, &outStream);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t20\tPARALOGY1REV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t60\t80\tREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "Sequence\t0\t20\tPARALOGY2\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 3);
*/
}

void WiggleLiftoverTest::testMultiBranchLifts(AlignmentConstPtr alignment)
{
}

void WiggleLiftoverTest::createCallBack(AlignmentPtr alignment)
{
  setupSharedAlignment(alignment);
}

void WiggleLiftoverTest::checkCallBack(AlignmentConstPtr alignment)
{
  testOneBranchLifts(alignment);
  testMultiBranchLifts(alignment);
}

void halBedLiftoverTest(CuTest *testCase)
{
  try
  {
    BedLiftoverTest tester;
    tester.check(testCase);
  }
  catch (...)
  {
    CuAssertTrue(testCase, false);
  }
}

void halWiggleLiftoverTest(CuTest *testCase)
{
  try
  {
    WiggleLiftoverTest tester;
    tester.check(testCase);
  }
  catch (...)
  {
    CuAssertTrue(testCase, false);
  }
}

CuSuite* halLiftoverTestSuite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halBedLiftoverTest);
  SUITE_ADD_TEST(suite, halWiggleLiftoverTest);
  return suite;
}

int halLiftoverRunAllTests(void) {
   CuString *output = CuStringNew();
   CuSuite* suite = CuSuiteNew();
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

