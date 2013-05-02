/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <ctime>
#include "halMappedSegmentTest.h"
#include "halRandomData.h"
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"
#include "hal.h"


using namespace std;
using namespace hal;

void MappedSegmentMapUpTest::createCallBack(AlignmentPtr alignment)
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

void MappedSegmentMapUpTest::testTopSegment(AlignmentConstPtr alignment,
                                            TopSegmentIteratorConstPtr top)
{
  const Genome* parent = top->getGenome()->getParent();
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, parent, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  CuAssertTrue(_testCase, mseg->getGenome() == bottom->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == bottom->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == bottom->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == bottom->getReversed());
}

void MappedSegmentMapUpTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");
  TopSegmentIteratorConstPtr top = child2->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
  top = child1->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
}

void MappedSegmentMapDownTest::testBottomSegment(
  AlignmentConstPtr alignment,
  BottomSegmentIteratorConstPtr bottom,
  hal_size_t childIndex)
{
  const Genome* child = bottom->getGenome()->getChild(childIndex);
  vector<MappedSegmentConstPtr> results;
  bottom->getMappedSegments(results, child, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == 
               bottom->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               bottom->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == bottom->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == bottom->getReversed());
  TopSegmentIteratorConstPtr top = child->getTopSegmentIterator();
  top->toChild(bottom, childIndex);
  CuAssertTrue(_testCase, mseg->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == top->getReversed());
}

void MappedSegmentMapDownTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");

  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
  bottom->slice(1,2);
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
  bottom->toReverse();
  testBottomSegment(alignment, bottom, 0);
  testBottomSegment(alignment, bottom, 1);
}

void MappedSegmentMapAcrossTest::testTopSegment(AlignmentConstPtr alignment,
                                                TopSegmentIteratorConstPtr top)
{
  const Genome* parent = top->getGenome()->getParent();
  const Genome* other = top->getGenome()->getName() == "child1" ? 
     alignment->openGenome("child2") : alignment->openGenome("child1");
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, other, NULL, false);
  CuAssertTrue(_testCase, results.size() == 1);
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  TopSegmentIteratorConstPtr sister = other->getTopSegmentIterator();
  sister->toChildG(bottom, other);
  CuAssertTrue(_testCase, mseg->getGenome() == sister->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == sister->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == sister->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == sister->getReversed());
}

void MappedSegmentMapAcrossTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");
  TopSegmentIteratorConstPtr top = child2->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
  top = child1->getTopSegmentIterator();
  testTopSegment(alignment, top);
  top->slice(1,2);
  testTopSegment(alignment, top);
  top->toReverse();
  testTopSegment(alignment, top);
}

void MappedSegmentMapDupeTest::createCallBack(AlignmentPtr alignment)
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
  seqVec[0] = Sequence::Info("Sequence", 3, 0, 1);
  parent->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  child1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 9, 3, 0);
  child2->setDimensions(seqVec);

  parent->setString("CCC");
  child1->setString("CCCTACGTG");
  child2->setString("CCCTACGTG");

  bi = parent->getBottomSegmentIterator();
  bs.set(0, 3, 0);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs._children.push_back(pair<hal_size_t, bool>(0, false));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 3, 0, true, 0, 1);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, true, 0, 2);
  ti->toRight();
  ts.set(6, 3, 0, true, 0, 0);

  ti = child2->getTopSegmentIterator();
  ts.set(0, 3, 0, true, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 3, 0, true, 0);
  ti->toRight();
  ts.set(6, 3, 0, true, 0);
}

void MappedSegmentMapDupeTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* parent = alignment->openGenome("parent");
  const Genome* child1 = alignment->openGenome("child1");
  const Genome* child2 = alignment->openGenome("child2");

  TopSegmentIteratorConstPtr top = child1->getTopSegmentIterator();
  vector<MappedSegmentConstPtr> results;
  top->getMappedSegments(results, child2, NULL, true);
  CuAssertTrue(_testCase, results.size() == 1);
  
  MappedSegmentConstPtr mseg = results[0];
  CuAssertTrue(_testCase, mseg->getSource()->getGenome() == top->getGenome());
  CuAssertTrue(_testCase, mseg->getSource()->getStartPosition() == 
               top->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getLength() == top->getLength());
  CuAssertTrue(_testCase, 
               mseg->getSource()->getReversed() == top->getReversed());
  BottomSegmentIteratorConstPtr bottom = parent->getBottomSegmentIterator();
  bottom->toParent(top);
  TopSegmentIteratorConstPtr sister = child2->getTopSegmentIterator();
  sister->toChildG(bottom, child2);
  CuAssertTrue(_testCase, mseg->getGenome() == sister->getGenome());
  CuAssertTrue(_testCase, 
               mseg->getStartPosition() == sister->getStartPosition());
  CuAssertTrue(_testCase, 
               mseg->getLength() == sister->getLength());
  CuAssertTrue(_testCase, 
               mseg->getReversed() == sister->getReversed());
}

void halMappedSegmentMapUpTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapUpTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapDownTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapDownTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapAcrossTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapAcrossTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halMappedSegmentMapDupeTest(CuTest *testCase)
{
  try 
  {
    MappedSegmentMapDupeTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}


CuSuite* halMappedSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halMappedSegmentMapUpTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDownTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapAcrossTest);
  SUITE_ADD_TEST(suite, halMappedSegmentMapDupeTest);
  return suite;
}

