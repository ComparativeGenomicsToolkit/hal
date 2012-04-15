/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <H5Cpp.h>
#include <H5Exception.h>
#include "halGenomeTest.h"
#include "halAlignmentTest.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignment.h"
#include "halGenome.h"
#include "halGenome.h"
#include "halMetaData.h"
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

void GenomeMetaTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
 
  MetaData* ancMeta = ancGenome->getMetaData();
  ancMeta->set("Young", "Jeezy");
}

void GenomeMetaTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openConstGenome("AncGenome");
  const MetaData* ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 0);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 0);
  CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == 0);
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
}

void GenomeCreateTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome* leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);
  Genome* leaf2Genome = alignment->addLeafGenome("Leaf2", "AncGenome", 0.2);
  Genome* leaf3Genome = alignment->addLeafGenome("Leaf3", "AncGenome", 0.3);
 
  MetaData* ancMeta = ancGenome->getMetaData();
  ancMeta->set("Young", "Jeezy");

  
  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 700000, 0);
  leaf1Genome->setDimensions(seqVec);
  seqVec[0] =Sequence::Info("Sequence", 2000000, 700000, 0);
  leaf2Genome->setDimensions(seqVec);
  seqVec[0] =Sequence::Info("Sequence", 3000000, 700000, 0);
  leaf3Genome->setDimensions(seqVec);
}

void GenomeCreateTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openConstGenome("AncGenome");
  const MetaData* ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  const Genome* leaf1Genome = alignment->openConstGenome("Leaf1");
  const Genome* leaf2Genome = alignment->openConstGenome("Leaf2");
  const Genome* leaf3Genome = alignment->openConstGenome("Leaf3");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, leaf1Genome->getName() == "Leaf1");
  CuAssertTrue(_testCase, leaf2Genome->getName() == "Leaf2");  
  CuAssertTrue(_testCase, leaf3Genome->getName() == "Leaf3");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 5000);
  CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == 700000);
  CuAssertTrue(_testCase, leaf1Genome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, leaf1Genome->getNumTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf1Genome->getNumBottomSegments() == 0);
  CuAssertTrue(_testCase, leaf2Genome->getSequenceLength() == 2000000);
  CuAssertTrue(_testCase, leaf2Genome->getNumTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf2Genome->getNumBottomSegments() == 0);
  CuAssertTrue(_testCase, leaf3Genome->getSequenceLength() == 3000000);
  CuAssertTrue(_testCase, leaf3Genome->getNumTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf3Genome->getNumBottomSegments() == 0);
}

void GenomeUpdateTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);  
  alignment->close();

  alignment->open(_createPath, false);
  ancGenome = alignment->openGenome("AncGenome");
  seqVec[0] = Sequence::Info("Sequence", 10000005, 14000, 2000001);
  ancGenome->setDimensions(seqVec);  
}

void GenomeUpdateTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openConstGenome("AncGenome");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 10000005);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 14000);
  CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == 2000001);
}

void halGenomeMetaTest(CuTest *testCase)
{
  GenomeMetaTest tester;
  tester.check(testCase);
}

void halGenomeCreateTest(CuTest *testCase)
{
  try
  {
    GenomeCreateTest tester;
    tester.check(testCase);
  }
  catch (H5::Exception& e) 
  {
    cerr << e.getDetailMsg() << " " << endl;
    CuAssertTrue(testCase, false);
  }
}

void halGenomeUpdateTest(CuTest *testCase)
{
  try
  {
    GenomeUpdateTest tester;
    tester.check(testCase);
  }
  catch (H5::Exception& e) 
  {
    cerr << e.getDetailMsg() << " " << endl;
    CuAssertTrue(testCase, false);
  }
}


CuSuite* halGenomeTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halGenomeMetaTest);
  SUITE_ADD_TEST(suite, halGenomeCreateTest);
  SUITE_ADD_TEST(suite, halGenomeUpdateTest);
  return suite;
}

