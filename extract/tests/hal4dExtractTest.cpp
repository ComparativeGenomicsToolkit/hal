#include "hal.h"
#include "hal4dExtract.h"
#include "hal4dExtractTest.h"

using namespace std;
using namespace hal;

void ConservedBed4dExtractTest::createCallBack(AlignmentPtr alignment)
{
  Genome *root = alignment->addRootGenome("root");
  Genome *leaf1 = alignment->addLeafGenome("leaf1", "root", 1);
  Genome *leaf2 = alignment->addLeafGenome("leaf2", "root", 1);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("rootSequence", 192, 0, 1);
  root->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("leaf1Sequence", 192, 1, 0);
  leaf1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("leaf2Sequence", 192, 1, 0);
  leaf2->setDimensions(seqVec);
  // all possible codons, in a single frame -- we will test only a few but make
  // sure that the correct number of sites is output as well
  // changed final 4d site, which should change nothing since ancestors are
  // ignored
  root->setString("aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgttttgatgctggtgtttattcttgttt");
  BottomSegmentIteratorPtr botIt = root->getBottomSegmentIterator();
  botIt->setCoordinates(0, 192);
  botIt->setChildIndex(0, 0);
  botIt->setChildIndex(1, 0);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);

  // added unconserved 4d site at 2, disabled 4d site at 14 and 18/20
  leaf1->setString("acaaacaagaataaaaccatgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
  TopSegmentIteratorPtr topIt = leaf1->getTopSegmentIterator();
  topIt->setCoordinates(0, 192);
  topIt->setParentIndex(0);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);

  // disabled 4d site at 14
  leaf2->setString("aaaaacaagaataaaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
  topIt = leaf2->getTopSegmentIterator();
  topIt->setCoordinates(0, 192);
  topIt->setParentIndex(0);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
}

void ConservedBed4dExtractTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome *genome = alignment->openGenome("root");
  stringstream bedFile("rootSequence\t0\t192\tFORWARD\t0\t+\n"
                       "rootSequence\t0\t192\tREV\t0\t-\n");
  stringstream outStream;
  Extract4d extract;
  extract.run(genome, &bedFile, &outStream, -1, true);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t14\t15\tFORWARD\t0\t+") == streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t177\t178\tREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t18\t19\tREV\t0\t-") == streamResults.end());
  // 14, 18(+ strand) /20 (- strand) disabled, so 3 are missing
  CuAssertTrue(_testCase, streamResults.size() == 61);
}

void Bed4dExtractTest::createCallBack(AlignmentPtr alignment)
{
  Genome *genome = alignment->addRootGenome("root");
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("rootSequence", 192, 0, 0);
  genome->setDimensions(seqVec);
  // all possible codons, in a single frame -- we will test only a few but make
  // sure that the correct number of sites is output as well
  genome->setString("aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
}

void Bed4dExtractTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome *genome = alignment->openGenome("root");
  stringstream bedFile("rootSequence\t0\t192\tFORWARD\t0\t+\n"
                       "rootSequence\t0\t192\tREV\t0\t-\n");
  stringstream outStream;
  Extract4d extract;
  extract.run(genome, &bedFile, &outStream, -1);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t14\t15\tFORWARD\t0\t+") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t177\t178\tREV\t0\t-") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 64);
}

void Block4dExtractTest::createCallBack(AlignmentPtr alignment)
{
  Genome *genome = alignment->addRootGenome("root");
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("rootSequence", 192, 0, 0);
  genome->setDimensions(seqVec);
  // all possible codons, in a single frame -- we will test only a few but make
  // sure that the correct number of sites is output as well
  genome->setString("aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
}

void Block4dExtractTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome *genome = alignment->openGenome("root");
  // test frame shift
  stringstream bedFile("rootSequence\t0\t192\tFORWARD\t0\t+\t0\t192\t0\t3\t17,6,7\t0,30,60\n"
                       "rootSequence\t0\t192\tREV\t0\t-\t0\t192\t0\t2\t13,17\t0,175");
  stringstream outStream;
  Extract4d extract;
  extract.run(genome, &bedFile, &outStream, -1);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t0\t192\tFORWARD\t0\t+\t0\t192\t0,0,0\t5\t1,1,1,1,1\t14,30,33,60,66") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t0\t192\tREV\t0\t-\t0\t192\t0,0,0\t4\t1,1,1,1\t3,6,12,177") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 2);
}

void ConservedBlock4dExtractTest::createCallBack(AlignmentPtr alignment)
{
  Genome *root = alignment->addRootGenome("root");
  Genome *leaf1 = alignment->addLeafGenome("leaf1", "root", 1);
  Genome *leaf2 = alignment->addLeafGenome("leaf2", "root", 1);
  vector<Sequence::Info> seqVec(2);
  seqVec[0] = Sequence::Info("rootSequence", 192, 0, 1);
  seqVec[1] = Sequence::Info("rootSequence2", 192, 0, 1);
  root->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("leaf1Sequence", 192, 1, 0);
  seqVec[1] = Sequence::Info("leaf1Sequence2", 192, 1, 0);
  leaf1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("leaf2Sequence", 192, 1, 0);
  seqVec[1] = Sequence::Info("leaf2Sequence2", 192, 1, 0);
  leaf2->setDimensions(seqVec);
  // all possible codons, in a single frame -- we will test only a few but make
  // sure that the correct number of sites is output as well
  // changed final 4d site, which should change nothing since ancestors are
  // ignored
  root->setString("aaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgttttgatgctggtgtttattcttgtttaaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgttttgatgctggtgtttattcttgttt");
  BottomSegmentIteratorPtr botIt = root->getBottomSegmentIterator();
  botIt->setCoordinates(0, 192);
  botIt->setChildIndex(0, 0);
  botIt->setChildIndex(1, 0);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);
  botIt->toRight();
  botIt->setCoordinates(192, 192);
  botIt->setChildIndex(0, 1);
  botIt->setChildIndex(1, 1);
  botIt->setChildReversed(0, false);
  botIt->setChildReversed(1, false);
  botIt->setTopParseIndex(NULL_INDEX);

  // added unconserved 4d site at 2, disabled 4d site at 14 and 18/20
  leaf1->setString("acaaacaagaataaaaccatgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgtttacaaacaagaataaaaccatgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
  TopSegmentIteratorPtr topIt = leaf1->getTopSegmentIterator();
  topIt->setCoordinates(0, 192);
  topIt->setParentIndex(0);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(192, 192);
  topIt->setParentIndex(1);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);

  // disabled 4d site at 14
  leaf2->setString("aaaaacaagaataaaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgtttaaaaacaagaataaaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
  topIt = leaf2->getTopSegmentIterator();
  topIt->setCoordinates(0, 192);
  topIt->setParentIndex(0);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
  topIt->toRight();
  topIt->setCoordinates(192, 192);
  topIt->setParentIndex(1);
  topIt->setParentReversed(false);
  topIt->setNextParalogyIndex(NULL_INDEX);
  topIt->setBottomParseIndex(NULL_INDEX);
}

void ConservedBlock4dExtractTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome *genome = alignment->openGenome("root");
  // test frame shift
  stringstream bedFile("rootSequence\t0\t192\tFORWARD\t0\t+\t0\t192\t0\t3\t17,6,7\t0,30,60\n"
                       "rootSequence\t0\t192\tREV\t0\t-\t0\t192\t0\t2\t13,17\t0,175\n"
                       "rootSequence2\t0\t192\tFORWARD\t0\t+\t0\t192\t0\t3\t17,6,7\t0,30,60\n"
                       "rootSequence2\t0\t192\tREV\t0\t-\t0\t192\t0\t2\t13,17\t0,175\n");
  stringstream outStream;
  Extract4d extract;
  extract.run(genome, &bedFile, &outStream, -1, true);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t0\t192\tFORWARD\t0\t+\t0\t192\t0,0,0\t2\t1,1\t33,66") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t0\t192\tREV\t0\t-\t0\t192\t0,0,0\t3\t1,1,1\t3,6,177") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence2\t0\t192\tFORWARD\t0\t+\t0\t192\t0,0,0\t2\t1,1\t33,66") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence2\t0\t192\tREV\t0\t-\t0\t192\t0,0,0\t3\t1,1,1\t3,6,177") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 4);
}

void CDS4dExtractTest::createCallBack(AlignmentPtr alignment)
{
  Genome *genome = alignment->addRootGenome("root");
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("rootSequence", 213, 0, 0);
  genome->setDimensions(seqVec);
  // all possible codons, in a single frame, starting at position 21
  // -- we will test only a few but make sure that the correct number
  // of sites is output as well
  genome->setString("caccatcatcatcatcatcataaaaacaagaatacaaccacgactagaagcaggagtataatcatgattcaacaccagcatccacccccgcctcgacgccggcgtctactcctgcttgaagacgaggatgcagccgcggctggaggcgggggtgtagtcgtggtttaatactagtattcatcctcgtcttgatgctggtgtttattcttgttt");
}

void CDS4dExtractTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome *genome = alignment->openGenome("root");
  stringstream bedFile("rootSequence\t1\t212\tFORWARD\t0\t+\t21\t88\t0\t5\t10,19,6,8,10\t0,18,50,80,90\n"
                       "rootSequence\t1\t213\tREV\t0\t-\t25\t198\t0\t2\t13,17\t20,195");
  stringstream outStream;
  Extract4d extract;
  extract.run(genome, &bedFile, &outStream, -1);
  vector<string> streamResults = chopString(outStream.str(), "\n");
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t1\t212\tFORWARD\t0\t+\t21\t88\t0,0,0\t5\t1,1,1,1,1\t34,50,53,80,86") != streamResults.end());
  CuAssertTrue(_testCase, find(streamResults.begin(), streamResults.end(), "rootSequence\t1\t213\tREV\t0\t-\t25\t198\t0,0,0\t2\t1,1\t26,32") != streamResults.end());
  CuAssertTrue(_testCase, streamResults.size() == 2);
}

void halBlock4dExtractTest(CuTest *testCase)
{
  // try
  // {
    Block4dExtractTest tester;
    tester.check(testCase);
  // }
  // catch (...)
  // {
  //   CuAssertTrue(testCase, false);
  // }
}

void halCDS4dExtractTest(CuTest *testCase)
{
  // try
  // {
    CDS4dExtractTest tester;
    tester.check(testCase);
  // }
  // catch (...)
  // {
  //   CuAssertTrue(testCase, false);
  // }
}

void halConservedBlock4dExtractTest(CuTest *testCase)
{
  // try
  // {
    ConservedBlock4dExtractTest tester;
    tester.check(testCase);
  // }
  // catch (...)
  // {
  //   CuAssertTrue(testCase, false);
  // }
}

CuSuite* hal4dExtractTestSuite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halBlock4dExtractTest);
  SUITE_ADD_TEST(suite, halConservedBlock4dExtractTest);
  SUITE_ADD_TEST(suite, halCDS4dExtractTest);
  return suite;
}

int hal4dExtractRunAllTests(void)
{
   CuString *output = CuStringNew();
   CuSuite* suite = CuSuiteNew();
   CuSuiteAddSuite(suite, hal4dExtractTestSuite());
   CuSuiteRun(suite);
   CuSuiteSummary(suite, output);
   CuSuiteDetails(suite, output);
   printf("%s\n", output->buffer);
   return suite->failCount > 0;
}

int main(int argc, char *argv[])
{
  return hal4dExtractRunAllTests();
}
