/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include "halGenomeTest.h"
#include "halAlignment.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignmentTest.h"
#include "halBottomSegmentIterator.h"
#include "halDNAIterator.h"
#include "halGenome.h"
#include "halGenome.h"
#include "halMetaData.h"
#include "halTopSegmentIterator.h"
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
  const Genome* ancGenome = alignment->openGenome("AncGenome");
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
  const Genome* dudGenome = alignment->openGenome("Zebra");
  CuAssertTrue(_testCase, dudGenome == NULL);
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  const MetaData* ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  const Genome* leaf1Genome = alignment->openGenome("Leaf1");
  const Genome* leaf2Genome = alignment->openGenome("Leaf2");
  const Genome* leaf3Genome = alignment->openGenome("Leaf3");
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
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 10000005);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 14000);
  CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == 2000001);
}

void GenomeStringTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  hal_size_t seqLength = 28889943;
  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", seqLength, 5000, 700000);
  ancGenome->setDimensions(seqVec);  
  
  _string = randomString(seqLength);
  ancGenome->setString(_string);
}

void GenomeStringTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  string genomeString;
  ancGenome->getString(genomeString);
  CuAssertTrue(_testCase, genomeString == _string);
}

void GenomeCopyTest::createCallBack(AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);

  // Hacky: Need a different alignment to test copying the bottom
  // segments correctly.  (the names of a node's children are used
  // when copying bottom segments, and two genomes can't have the same
  // name in the same alignment)
  vector<AlignmentPtr> createInstances = getTestAlignmentInstances();
  _secondAlignment = createInstances[0];
  _path = getTempFile();
  _secondAlignment->createNew(_path);

  Genome* ancGenome = alignment->addRootGenome("AncGenome", 0);
  Genome *leafGenome = alignment->addLeafGenome("LeafGenome1",
                                                "AncGenome", 0);
  // This genome will test copyDimensions, copyTopSegments,
  // copyBottomSegments, copySequence, copyMetadata
  Genome* copyRootGenome = _secondAlignment->addRootGenome("copyRootGenome",
                                                           0);
  Genome *copyLeafGenome = _secondAlignment->addLeafGenome("LeafGenome1",
                                                           "copyRootGenome",
                                                           0);

  MetaData* ancMeta = ancGenome->getMetaData();
  ancMeta->set("Young", "Jeezy");

  vector<Sequence::Info> seqVec(1);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 0, 700000);
  ancGenome->setDimensions(seqVec);
  seqVec[0] =Sequence::Info("Sequence", 1000000, 5000, 0);
  leafGenome->setDimensions(seqVec);
  string ancSeq = "CAT";
  hal_index_t n = ancGenome->getSequenceLength();
  DNAIteratorPtr dnaIt = ancGenome->getDNAIterator();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % ancSeq.size();
    dnaIt->setChar(ancSeq[i]);
  }
  n = leafGenome->getSequenceLength();
  dnaIt = leafGenome->getDNAIterator();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % ancSeq.size();
    dnaIt->setChar(ancSeq[i]);
  }

  TopSegmentIteratorPtr topIt = leafGenome->getTopSegmentIterator();
  n = leafGenome->getNumTopSegments();
  for (; topIt->getArrayIndex() < n; topIt->toRight())
  {
    topIt->setCoordinates(topIt->getArrayIndex(), 1);
    topIt->setParentIndex(3);
    topIt->setParentReversed(true);
    topIt->setBottomParseIndex(5);
    if (topIt->getArrayIndex() != 6) {
      topIt->setNextParalogyIndex(6);
    } else {
      topIt->setNextParalogyIndex(7);
    }
  }
  BottomSegmentIteratorPtr botIt = ancGenome->getBottomSegmentIterator();
  n = ancGenome->getNumBottomSegments();
  for (; botIt->getArrayIndex() < n; botIt->toRight())
  {
    botIt->setCoordinates(botIt->getArrayIndex(), 1);
    botIt->setChildIndex(0, 3);
    botIt->setChildReversed(0, true);
    botIt->setTopParseIndex(5);
  }

  seqVec[0] =Sequence::Info("Sequence", 3300, 0, 1100);
  copyRootGenome->setDimensions(seqVec);
  seqVec[0] =Sequence::Info("Sequence", 3300, 2200, 0);
  copyLeafGenome->setDimensions(seqVec);
  string copySeq = "TAG";
  dnaIt = copyRootGenome->getDNAIterator();
  n = copyRootGenome->getSequenceLength();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % copySeq.size();
    dnaIt->setChar(copySeq[i]);
  }
  dnaIt = copyLeafGenome->getDNAIterator();
  n = copyLeafGenome->getSequenceLength();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % copySeq.size();
    dnaIt->setChar(copySeq[i]);
  }
  topIt = copyLeafGenome->getTopSegmentIterator();
  n = copyLeafGenome->getNumTopSegments();
  for (; topIt->getArrayIndex() < n; topIt->toRight())
  {
    topIt->setCoordinates(7, 8);
    topIt->setParentIndex(9);
    topIt->setParentReversed(false);
    topIt->setBottomParseIndex(11);
    if (topIt->getArrayIndex() != 12) {
      topIt->setNextParalogyIndex(12);
    } else {
      topIt->setNextParalogyIndex(7);
    }
  }
  botIt = copyRootGenome->getBottomSegmentIterator();
  n = copyRootGenome->getNumBottomSegments();
  for (; botIt->getArrayIndex() < n; botIt->toRight())
  {
    botIt->setCoordinates(6, 7);
    botIt->setChildIndex(0, 8);
    botIt->setChildReversed(0, false);
    botIt->setTopParseIndex(10);
  }
  
  ancGenome->copy(copyRootGenome);
  leafGenome->copy(copyLeafGenome);
  _secondAlignment->close();
}

void GenomeCopyTest::checkCallBack(hal::AlignmentConstPtr alignment)
{
  // FIXME: halAlignment->open() fails miserably but
  // openHalAlignmentReadOnly works? Probably some state isn't cleared
  // on close.
  AlignmentPtr tmp = hdf5AlignmentInstance();
  tmp->open(_path, true);
  _secondAlignment = tmp;
  const Genome* ancGenome = alignment->openGenome("AncGenome");
  CuAssertTrue(_testCase, ancGenome->getName() == "AncGenome");
  CuAssertTrue(_testCase, ancGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == 0);
  CuAssertTrue(_testCase, ancGenome->getNumBottomSegments() == 700000);
  const MetaData* ancMeta = ancGenome->getMetaData();
  CuAssertTrue(_testCase, ancMeta->get("Young") == "Jeezy");
  const Genome *leafGenome = alignment->openGenome("LeafGenome1");
  string ancSeq = "CAT";
  hal_index_t n = ancGenome->getSequenceLength();
  DNAIteratorConstPtr dnaIt = ancGenome->getDNAIterator();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % ancSeq.size();
    CuAssertTrue(_testCase, dnaIt->getChar() == ancSeq[i]);
  }
  TopSegmentIteratorPtr topIt = leafGenome->getTopSegmentIterator();
  n = leafGenome->getNumTopSegments();
  for (; topIt->getArrayIndex() < n; topIt->toRight())
  {
    CuAssertTrue(_testCase,
                 topIt->getStartPosition() == topIt->getArrayIndex());
    CuAssertTrue(_testCase, topIt->getLength() == 1);
    CuAssertTrue(_testCase, topIt->getParentIndex() == 3);
    CuAssertTrue(_testCase, topIt->getParentReversed() == true);
    CuAssertTrue(_testCase, topIt->getBottomParseIndex() == 5);
    if (topIt->getArrayIndex() != 6) {
      CuAssertTrue(_testCase, topIt->getNextParalogyIndex() == 6);
    } else {
      CuAssertTrue(_testCase, topIt->getNextParalogyIndex() == 7);
    }
  }
  BottomSegmentIteratorPtr botIt = ancGenome->getBottomSegmentIterator();
  n = ancGenome->getNumBottomSegments();
  for (; botIt->getArrayIndex() < n; botIt->toRight())
  {
    CuAssertTrue(_testCase,
                 botIt->getStartPosition() == botIt->getArrayIndex());
    CuAssertTrue(_testCase, botIt->getLength() == 1);
    CuAssertTrue(_testCase, botIt->getChildIndex(0) == 3);
    CuAssertTrue(_testCase, botIt->getChildReversed(0) == true);
    CuAssertTrue(_testCase, botIt->getTopParseIndex() == 5);
  }

  const Genome* copyRootGenome = _secondAlignment->openGenome("copyRootGenome");
  const Genome *copyLeafGenome = _secondAlignment->openGenome("LeafGenome1");
  CuAssertTrue(_testCase, copyRootGenome->getName() == "copyRootGenome");
  CuAssertTrue(_testCase, copyRootGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, copyRootGenome->getNumTopSegments() == 0);
  CuAssertTrue(_testCase, copyRootGenome->getNumBottomSegments() == 700000);
  CuAssertTrue(_testCase, copyLeafGenome->getName() == "LeafGenome1");
  CuAssertTrue(_testCase, copyLeafGenome->getSequenceLength() == 1000000);
  CuAssertTrue(_testCase, copyLeafGenome->getNumTopSegments() == 5000);
  CuAssertTrue(_testCase, copyLeafGenome->getNumBottomSegments() == 0);
  const MetaData* copyMeta = copyRootGenome->getMetaData();
  CuAssertTrue(_testCase, copyMeta->get("Young") == "Jeezy");
  n = copyRootGenome->getSequenceLength();
  dnaIt = copyRootGenome->getDNAIterator();
  for (; dnaIt->getArrayIndex() < n; dnaIt->toRight()) {
    size_t i = dnaIt->getArrayIndex() % ancSeq.size();
    CuAssertTrue(_testCase, dnaIt->getChar() == ancSeq[i]);
  }
  topIt = copyLeafGenome->getTopSegmentIterator();
  n = copyLeafGenome->getNumTopSegments();
  for (; topIt->getArrayIndex() < n; topIt->toRight())
  {
    CuAssertTrue(_testCase,
                 topIt->getStartPosition() == topIt->getArrayIndex());
    CuAssertTrue(_testCase, topIt->getLength() == 1);
    CuAssertTrue(_testCase, topIt->getParentIndex() == 3);
    CuAssertTrue(_testCase, topIt->getParentReversed() == true);
    CuAssertTrue(_testCase, topIt->getBottomParseIndex() == 5);
    if (topIt->getArrayIndex() != 6) {
      CuAssertTrue(_testCase, topIt->getNextParalogyIndex() == 6);
    } else {
      CuAssertTrue(_testCase, topIt->getNextParalogyIndex() == 7);
    }
  }
  botIt = copyRootGenome->getBottomSegmentIterator();
  n = copyRootGenome->getNumBottomSegments();
  for (; botIt->getArrayIndex() < n; botIt->toRight())
  {
    CuAssertTrue(_testCase,
                 botIt->getStartPosition() == botIt->getArrayIndex());
    CuAssertTrue(_testCase, botIt->getLength() == 1);
    CuAssertTrue(_testCase, botIt->getChildIndex(0) == 3);
    CuAssertTrue(_testCase, botIt->getChildReversed(0) == true);
    CuAssertTrue(_testCase, botIt->getTopParseIndex() == 5);
  }

  _secondAlignment->close();
  removeTempFile((char *)_path.c_str());
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
  catch (...) 
  {
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
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

void halGenomeStringTest(CuTest *testCase)
{
  try
  {
    GenomeStringTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

void halGenomeCopyTest(CuTest *testCase)
{
/*  try
    {*/
    GenomeCopyTest tester;
    tester.check(testCase);
/*  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
    } */
}

CuSuite* halGenomeTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halGenomeMetaTest);
  SUITE_ADD_TEST(suite, halGenomeCreateTest);
  SUITE_ADD_TEST(suite, halGenomeUpdateTest);
  SUITE_ADD_TEST(suite, halGenomeStringTest);
  SUITE_ADD_TEST(suite, halGenomeCopyTest);
  return suite;
}

