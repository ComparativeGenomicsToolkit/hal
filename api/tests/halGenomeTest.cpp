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

void GenomeMetaTest::createCallBack(hal::AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  GenomePtr ancGenome = alignment->addRootGenome("AncGenome", 0);
  GenomePtr leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);
  GenomePtr leaf2Genome = alignment->addLeafGenome("Leaf2", "AncGenome", 0.2);
  GenomePtr leaf3Genome = alignment->addLeafGenome("Leaf3", "AncGenome", 0.3);
 
  MetaDataPtr ancMeta = ancGenome->getMetaData();
  ancMeta->set("Young", "Jeezy");
}

void GenomeMetaTest::checkCallBack(hal::AlignmentConstPtr alignment)
{
  GenomeConstPtr ancGenome = alignment->openConstGenome("AncGenome");
  MetaDataConstPtr ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 0);
  CuAssertTrue(_testCase, ancGenome->getNumberTopSegments() == 0);
  CuAssertTrue(_testCase, ancGenome->getNumberBottomSegments() == 0);
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
}

void GenomeCreateTest::createCallBack(hal::AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  GenomePtr ancGenome = alignment->addRootGenome("AncGenome", 0);
  GenomePtr leaf1Genome = alignment->addLeafGenome("Leaf1", "AncGenome", 0.1);
  GenomePtr leaf2Genome = alignment->addLeafGenome("Leaf2", "AncGenome", 0.2);
  GenomePtr leaf3Genome = alignment->addLeafGenome("Leaf3", "AncGenome", 0.3);
 
  MetaDataPtr ancMeta = ancGenome->getMetaData();
  ancMeta->set("Young", "Jeezy");

  ancGenome->reset(1000000, 5000, 700000);
  leaf1Genome->reset(1000000, 700000, 0);
  leaf2Genome->reset(2000000, 700000, 0);
  leaf3Genome->reset(3000000, 700000, 0);
}

void GenomeCreateTest::checkCallBack(hal::AlignmentConstPtr alignment)
{
  GenomeConstPtr ancGenome = alignment->openConstGenome("AncGenome");
  MetaDataConstPtr ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  GenomeConstPtr leaf1Genome = alignment->openConstGenome("Leaf1");
  GenomeConstPtr leaf2Genome = alignment->openConstGenome("Leaf2");
  GenomeConstPtr leaf3Genome = alignment->openConstGenome("Leaf3");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, leaf1Genome->getName() == "Leaf1");
  CuAssertTrue(_testCase, leaf2Genome->getName() == "Leaf2");  
  CuAssertTrue(_testCase, leaf3Genome->getName() == "Leaf3");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, ancGenome->getNumberTopSegments() == 5000);
  CuAssertTrue(_testCase, ancGenome->getNumberBottomSegments() == 700000);
  CuAssertTrue(_testCase, leaf1Genome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, leaf1Genome->getNumberTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf1Genome->getNumberBottomSegments() == 0);
  CuAssertTrue(_testCase, leaf2Genome->getSequenceLength() == 2000000);
  CuAssertTrue(_testCase, leaf2Genome->getNumberTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf2Genome->getNumberBottomSegments() == 0);
  CuAssertTrue(_testCase, leaf3Genome->getSequenceLength() == 3000000);
  CuAssertTrue(_testCase, leaf3Genome->getNumberTopSegments() == 700000);
  CuAssertTrue(_testCase, leaf3Genome->getNumberBottomSegments() == 0);
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

CuSuite* halGenomeTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halGenomeMetaTest);
  SUITE_ADD_TEST(suite, halGenomeCreateTest);
  return suite;
}

