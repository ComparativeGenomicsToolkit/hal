/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <sstream>
#include "halGenomeTest.h"
#include "halAlignmentTest.h"
#include "halSequenceTest.h"

extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;


void SequenceCreateTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  
  size_t numSequences = 1000;
  vector<Sequence::Info> seqVec;
  for (size_t i = 0; i < numSequences; ++i)
  {
    std::stringstream ss;
    ss << i;
    hal_size_t len = 1 + i * 5 + i;
    string name = "sequence" + ss.str();
    seqVec.push_back(Sequence::Info(name, len, i, i * 2));
  }
  ancGenome->setDimensions(seqVec);
}

void SequenceCreateTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  
  hal_size_t numSequences = ancGenome->getNumSequences();
  CuAssertTrue(_testCase, numSequences = 1000);

  SequenceIteratorConstPtr seqIt = ancGenome->getSequenceIterator();
  for (hal_size_t i = 0; i < numSequences; ++i)
  {
    stringstream ss;
    ss << i;
    hal_size_t len = 1 + i * 5 + i;
    string name = "sequence" + ss.str();
    const Sequence* seq = seqIt->getSequence();
    CuAssertTrue(_testCase, seq->getName() == name);
    CuAssertTrue(_testCase, seq->getSequenceLength() == len);
    CuAssertTrue(_testCase, seq->getNumTopSegments() == i);
    CuAssertTrue(_testCase, seq->getNumBottomSegments() == i * 2);
    const Genome* gen = seq->getGenome();
    CuAssertTrue(_testCase, gen->getName() == "AncGenome");
    seqIt->toNext();
  }

  const Sequence* seq = ancGenome->getSequence("sequence555");
  CuAssertTrue(_testCase, seq->getName() == "sequence555");
  seq = ancGenome->getSequenceBySite(0);
  CuAssertTrue(_testCase, seq->getName() == "sequence0");
  seq = ancGenome->getSequenceBySite(45);
  CuAssertTrue(_testCase, seq->getName() == "sequence4");
}


void halSequenceCreateTest(CuTest *testCase)
{
  try
  {
    SequenceCreateTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}


CuSuite* halSequenceTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halSequenceCreateTest);
  return suite;
}

