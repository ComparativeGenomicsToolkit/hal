/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halChainTests.h"
#include "halBlockViz.h"
#include "halBlockMapper.h"
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"
#include "halChainGetBlocksTest.h"

using namespace std;
using namespace hal;


void ChainGetBlocksSimpleTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // setup simple case were there is an edge from a parent to 
  // child and it is reversed
  Genome* parent = alignment->addRootGenome("parent");
  Genome* child1 = alignment->addLeafGenome("child1", "parent", 1);
  Genome* child2 = alignment->addLeafGenome("child2", "parent", 1);
  seqVec[0] = Sequence::Info("Sequence", 10, 0, 1);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 10, 1, 0);
  child1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 10, 1, 0);
  child2->setDimensions(seqVec);

  parent->setString("CCCTACGTGC");
  child1->setString("CCCTACGTGC");
  child2->setString("CCCTACGTGC");

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 10, 0);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 10, 0, true, 0);
  ts.applyTo(ti);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 10, 0, false, 0);
  ts.applyTo(ti);
}

void ChainGetBlocksSimpleTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child2");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(parent, child, 0, 9, false, false, 0, true);
  mapper.map();
  const BlockMapper::MSSet& segMap = mapper.getMap();
  CuAssertTrue(_testCase, segMap.size() == 1);
  MSRefSet::const_iterator mapIt = segMap.begin();
  
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 10);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 10);    
}

void ChainGetBlocksInversionTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child1");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(parent, child, 0, 9, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 1);
  MSRefSet::const_iterator mapIt = segMap.begin();
  
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 10);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 10);    
}

void ChainGetBlocksOffsetTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child2");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(parent, child, 1, 4, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 3);
  
  MSRefSet::const_iterator mapIt = segMap.begin();
  
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 1);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 1);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 4);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 4);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 5);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 5);    

}

void ChainGetBlocksInversionOffsetTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child1");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(parent, child, 1, 4, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 3);
  
  MSRefSet::const_iterator mapIt = segMap.begin();
 
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 1);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 1);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 4);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 8);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 4);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 5);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 4);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 5);    

}

void ChainGetBlocksOffsetQRefTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child2");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(child, parent, 1, 4, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 3);
  
  MSRefSet::const_iterator mapIt = segMap.begin();
 
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 1);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 1);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 4);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 4);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 5);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, queSeg->getReversed() == false);  
  CuAssertTrue(_testCase, queSeg->getLength() == 5);    

}

void ChainGetBlocksInversionOffsetQRefTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child1");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(child, parent, 1, 4, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 3);
  
  MSRefSet::const_iterator mapIt = segMap.begin();
 
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 1);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 1);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 4);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 8);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 4);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 5);  

  CuAssertTrue(_testCase, queSeg->getGenome() == parent);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 4);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 5);    

}

 
void ChainGetBlocksInversionOffsetQSisTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(child1, child2, 1, 4, false, false, 0, true);
  mapper.map();
  MSRefSet segMap;
  segMap.insert(mapper.getMap().begin(), mapper.getMap().end());
  CuAssertTrue(_testCase, segMap.size() == 3);
  
  MSRefSet::const_iterator mapIt = segMap.begin();
 
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child1);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 0);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 1);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child2);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 1);    

  ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child1);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 1);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 4);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child2);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 8);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 4);    

 ++mapIt;
  refSeg = (*mapIt)->getSource();
  queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == child1);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 5);  
  CuAssertTrue(_testCase, refSeg->getReversed() == false);  
  CuAssertTrue(_testCase, refSeg->getLength() == 5);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child2);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 4);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 5);    

}

void ChainGetBlocksSimpleLiftoverTest::checkCallBack(
  AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child = alignment->openGenome("child2");

  // BASIC TEST
  BlockMapper mapper;
  mapper.init(parent, child, 0, 9, true, false, 0, false);
  mapper.map();
  const BlockMapper::MSSet& segMap = mapper.getMap();
  CuAssertTrue(_testCase, segMap.size() == 1);
  MSRefSet::const_iterator mapIt = segMap.begin();
  
  SlicedSegmentConstPtr refSeg = (*mapIt)->getSource();
  SlicedSegmentConstPtr queSeg = (*mapIt);

  CuAssertTrue(_testCase, refSeg->getGenome() == parent);  
  CuAssertTrue(_testCase, refSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, refSeg->getReversed() == true);  
  CuAssertTrue(_testCase, refSeg->getLength() == 10);  

  CuAssertTrue(_testCase, queSeg->getGenome() == child);
  CuAssertTrue(_testCase, queSeg->getStartPosition() == 9);  
  CuAssertTrue(_testCase, queSeg->getReversed() == true);  
  CuAssertTrue(_testCase, queSeg->getLength() == 10);    
}

void halChainGetBlocksSimpleTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksSimpleTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksInversionTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksInversionTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksOffsetTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksOffsetTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksInversionOffsetTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksInversionOffsetTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksOffsetQRefTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksOffsetQRefTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksInversionOffsetQRefTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksInversionOffsetQRefTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksInversionOffsetQSisTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksInversionOffsetQSisTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halChainGetBlocksSimpleLiftoverTest(CuTest *testCase)
{
  try 
  {
    ChainGetBlocksSimpleLiftoverTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}


CuSuite *halChainGetBlocksTestSuite(void)
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halChainGetBlocksSimpleTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksInversionTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksOffsetTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksInversionOffsetTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksOffsetQRefTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksInversionOffsetQRefTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksInversionOffsetQSisTest);
  SUITE_ADD_TEST(suite, halChainGetBlocksSimpleLiftoverTest);
  return suite;
}

