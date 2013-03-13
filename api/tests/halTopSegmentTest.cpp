/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "halTopSegmentTest.h"
#include "halBottomSegmentTest.h"

using namespace std;
using namespace hal;

// just create a bunch of garbage data.  we don't care 
// about logical consistency for this test, just whether or not
// it's read and written properly. 
void TopSegmentStruct::setRandom()
{
  _length = rand();
  _startPosition = rand();
  _nextParalogyIndex = rand();
  _parentIndex = rand();
  _arrayIndex = rand();
  _bottomParseIndex = rand();
}

void TopSegmentStruct::set(hal_index_t startPosition,
                           hal_size_t length,
                           hal_index_t parentIndex,
                           bool parentReversed,
                           hal_index_t bottomParseIndex,
                           hal_index_t nextParalogyIndex)
{
  _startPosition = startPosition;
  _length = length;
  _parentIndex = parentIndex;
  _parentReversed = parentReversed;
  _bottomParseIndex = bottomParseIndex;
  _nextParalogyIndex = nextParalogyIndex;
}

void TopSegmentStruct::applyTo(TopSegmentIteratorPtr it) const
{
  TopSegment* seg = it->getTopSegment();
  seg->setCoordinates(_startPosition, _length);
  seg->setNextParalogyIndex(_nextParalogyIndex);
  seg->setParentIndex(_parentIndex);
  seg->setParentReversed(_parentReversed);
  seg->setBottomParseIndex(_bottomParseIndex);
}

void TopSegmentStruct::compareTo(TopSegmentIteratorConstPtr it, 
                                 CuTest* testCase) const
{
  const TopSegment* seg = it->getTopSegment();
  CuAssertTrue(testCase, _length == seg->getLength());
  CuAssertTrue(testCase, _startPosition == seg->getStartPosition());
  CuAssertTrue(testCase, _nextParalogyIndex == seg->getNextParalogyIndex());
  CuAssertTrue(testCase, _parentIndex == seg->getParentIndex());
  CuAssertTrue(testCase, _bottomParseIndex == seg->getBottomParseIndex());
}

void TopSegmentSimpleIteratorTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("Anc0", 0);
  size_t numChildren = 9;
  for (size_t i = 0; i < numChildren; ++i)
  {
    std::stringstream ss;
    ss << i;
    alignment->addLeafGenome(string("Leaf") + ss.str(), "Anc0", 0.1);
  }
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 10000);
  ancGenome->setDimensions(seqVec);
  
  CuAssertTrue(_testCase, ancGenome->getNumChildren() == numChildren);
  
  _topSegments.clear();
  for (size_t i = 0; i < ancGenome->getNumTopSegments(); ++i)
  {
    TopSegmentStruct topSeg;
    topSeg.setRandom();
    topSeg._length = 
       ancGenome->getSequenceLength() / ancGenome->getNumTopSegments();
    topSeg._startPosition = i * topSeg._length;
    _topSegments.push_back(topSeg);
  }
  
  TopSegmentIteratorPtr tsIt = ancGenome->getTopSegmentIterator(0);
  TopSegmentIteratorConstPtr tsEnd = 
     ancGenome->getTopSegmentEndIterator();
  for (size_t i = 0; tsIt != tsEnd; tsIt->toRight(), ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    _topSegments[i].applyTo(tsIt);
  }
}

void TopSegmentSimpleIteratorTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  CuAssertTrue(_testCase, 
               ancGenome->getNumTopSegments() == _topSegments.size());
  TopSegmentIteratorConstPtr tsIt = ancGenome->getTopSegmentIterator(0);
  for (size_t i = 0; i < ancGenome->getNumTopSegments(); ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    _topSegments[i].compareTo(tsIt, _testCase);
    tsIt->toRight();
  }
  tsIt = ancGenome->getTopSegmentIterator(
    ancGenome->getNumTopSegments() - 1);
  for (hal_index_t i = ancGenome->getNumTopSegments() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    _topSegments[i].compareTo(tsIt, _testCase);
    tsIt->toLeft();
  }

  tsIt = ancGenome->getTopSegmentIterator(0); 
  tsIt->slice(0, tsIt->getLength() - 1);
  for (hal_index_t i = 0; i < (hal_index_t)ancGenome->getSequenceLength(); ++i)
  {
    CuAssertTrue(_testCase, tsIt->getLength() == 1);
    CuAssertTrue(_testCase, tsIt->getStartPosition() == i);
    tsIt->toRight(tsIt->getStartPosition() + 1);
  }
  tsIt = ancGenome->getTopSegmentIterator(
    ancGenome->getNumTopSegments() - 1);
  tsIt->slice(tsIt->getLength() - 1, 0);
  for (hal_index_t i = ancGenome->getSequenceLength() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, tsIt->getLength() == 1);
    CuAssertTrue(_testCase, tsIt->getStartPosition() == i);
    tsIt->toLeft(tsIt->getStartPosition() - 1);
  }

  tsIt = ancGenome->getTopSegmentIterator(0); 
  tsIt->toReverse();
  CuAssertTrue(_testCase, tsIt->getReversed() == true);
  tsIt->slice(tsIt->getLength() - 1, 0);
  for (hal_index_t i = 0; i < (hal_index_t)ancGenome->getSequenceLength(); ++i)
  {
    CuAssertTrue(_testCase, tsIt->getLength() == 1);
    CuAssertTrue(_testCase, tsIt->getStartPosition() == i);
    tsIt->toLeft(tsIt->getStartPosition() + 1);
  }
  tsIt = ancGenome->getTopSegmentIterator(
    ancGenome->getNumTopSegments() - 1);
  tsIt->toReverse();
  tsIt->slice(0, tsIt->getLength() - 1);
  for (hal_index_t i = ancGenome->getSequenceLength() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, tsIt->getLength() == 1);
    CuAssertTrue(_testCase, tsIt->getStartPosition() == i);
    tsIt->toRight(tsIt->getStartPosition() - 1);
  }
}

void TopSegmentSequenceTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("Anc0", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);

  ancGenome->setSubString("CACACATTC", 500, 9);
  TopSegmentIteratorPtr tsIt = ancGenome->getTopSegmentIterator(100);
  tsIt->getTopSegment()->setCoordinates(500, 9);
}

void TopSegmentSequenceTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  TopSegmentIteratorConstPtr tsIt = ancGenome->getTopSegmentIterator(100);
  CuAssertTrue(_testCase, tsIt->getTopSegment()->getStartPosition() == 500);
  CuAssertTrue(_testCase, tsIt->getTopSegment()->getLength() == 9);
  string seq;
  tsIt->getString(seq);
  CuAssertTrue(_testCase, seq == "CACACATTC");
  tsIt->toReverse();
  tsIt->getString(seq);
  CuAssertTrue(_testCase, seq == "GAATGTGTG");
}

void TopSegmentIteratorParseTest::createCallBack(AlignmentPtr alignment)
{
 vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // case 1: bottom segment aligns perfectly with top segment
  Genome* case1 = alignment->addRootGenome("case1");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 2);
  case1->setDimensions(seqVec);
  
  ti = case1->getTopSegmentIterator();
  ts.set(0, 10, NULL_INDEX, false, 0, NULL_INDEX);
  ts.applyTo(ti);
  
  bi = case1->getBottomSegmentIterator();
  bs.set(0, 10, 0);
  bs.applyTo(bi);

  // case 2: bottom segment is completely contained in top segment
  Genome* case2 = alignment->addRootGenome("case2");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 3);
  case2->setDimensions(seqVec);
  
  ti = case2->getTopSegmentIterator();
  ts.set(0, 9, NULL_INDEX, false, 0, NULL_INDEX);
  ts.applyTo(ti);

  bi = case2->getBottomSegmentIterator();
  bs.set(0, 3, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(3, 4, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(7, 3, 0);
  bs.applyTo(bi);

  // case 3 top segment is completely contained in bottom segment
  Genome* case3 = alignment->addRootGenome("case3");
  seqVec[0] = Sequence::Info("Sequence", 10, 3, 2);
  case3->setDimensions(seqVec);

  ti = case3->getTopSegmentIterator();
  ts.set(0, 3, NULL_INDEX, false, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 4, NULL_INDEX, false, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(7, 3, NULL_INDEX, false, 0);
  ts.applyTo(ti);

  bi = case3->getBottomSegmentIterator();
  bs.set(0, 9, 0);
  bs.applyTo(bi);
 
  // case 4: top segment overhangs bottom segment on the left
  Genome* case4 = alignment->addRootGenome("case4");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 2);
  case4->setDimensions(seqVec);

  ti = case4->getTopSegmentIterator();
  ts.set(0, 9, NULL_INDEX, false, 0);
  ts.applyTo(ti);

  bi = case4->getBottomSegmentIterator();
  bs.set(0, 5, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(5, 5, 0);
  bs.applyTo(bi);
}

void TopSegmentIteratorParseTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  // case 1
  const Genome* case1 = alignment->openGenome("case1");
  ti = case1->getTopSegmentIterator();
  bi = case1->getBottomSegmentIterator();
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  CuAssertTrue(_testCase, bi->getLength() == ti->getLength());
  bi->slice(3, 1);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getLength() == 
               bi->getBottomSegment()->getLength() - 4);

  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  CuAssertTrue(_testCase, bi->getLength() == ti->getLength());

  // case 2
  const Genome* case2 = alignment->openGenome("case2");
  ti = case2->getTopSegmentIterator();
  bi = case2->getBottomSegmentIterator(1);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  bi->slice(1, 1);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());

  // case 3
  const Genome* case3 = alignment->openGenome("case3");
  ti = case3->getTopSegmentIterator();
  bi = case3->getBottomSegmentIterator();
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  bi->slice(2, 1);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());

  // case 4
  const Genome* case4 = alignment->openGenome("case4");
  ti = case4->getTopSegmentIterator();
  bi = case4->getBottomSegmentIterator(1);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  bi->slice(2, 2);
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
}

void TopSegmentIteratorToSiteTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // case 1: single segment
  Genome* case1 = alignment->addRootGenome("case1");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 0);
  case1->setDimensions(seqVec);
  ti = case1->getTopSegmentIterator();
  ts.set(0, 9);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(9, 1);
  ts.applyTo(ti);
  case1 = NULL;

  // case 2: bunch of random segments
  const hal_size_t numSegs = 1133;
  hal_size_t total = 0;
  vector<hal_size_t> segLens(numSegs);
  for (size_t i = 0 ; i < numSegs; ++i)
  {
    hal_size_t len = rand() % 77 + 1;
    segLens[i] = len;
    total += len;
    assert(len > 0);
  }
  Genome* case2 = alignment->addRootGenome("case2");
  seqVec[0] = Sequence::Info("Sequence", total, numSegs, 0);
  case2->setDimensions(seqVec);
  hal_index_t prev = 0;
  for (size_t i = 0 ; i < numSegs; ++i)
  {
    ti = case2->getTopSegmentIterator((hal_index_t)i);
    ts.set(prev, segLens[i]);
    prev += segLens[i];
    ts.applyTo(ti);
  }
}

void TopSegmentIteratorToSiteTest::checkGenome(const Genome* genome)
{
  TopSegmentIteratorConstPtr ti = genome->getTopSegmentIterator();
  for (hal_index_t pos = 0; 
       pos < (hal_index_t)genome->getSequenceLength(); ++pos)
  {
    ti->toSite(pos);
    CuAssertTrue(_testCase, ti->getStartPosition() == pos);
    CuAssertTrue(_testCase, ti->getLength() == 1);
    ti->toSite(pos, false);
    CuAssertTrue(_testCase, pos >= ti->getStartPosition() && 
                 pos < ti->getStartPosition() + (hal_index_t)ti->getLength());
    CuAssertTrue(_testCase, 
                 ti->getLength() == ti->getTopSegment()->getLength());
  }
}

void TopSegmentIteratorToSiteTest::checkCallBack(AlignmentConstPtr alignment)
{
  TopSegmentIteratorConstPtr bi;

  // case 1
  const Genome* case1 = alignment->openGenome("case1");
  checkGenome(case1);

  // case 2
  const Genome* case2 = alignment->openGenome("case2");
  checkGenome(case2);
}

void TopSegmentIteratorReverseTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  // setup simple case were there is an edge from a parent to 
  // child and it is reversed
  Genome* parent1 = alignment->addRootGenome("parent1");
  Genome* child1 = alignment->addLeafGenome("child1", "parent1", 1);
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 2);
  parent1->setDimensions(seqVec);
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 2);
  child1->setDimensions(seqVec);

  parent1->setString("CCCTACGTGC");
  child1->setString("CCCTACGTGC");

  bi = parent1->getBottomSegmentIterator();
  bs.set(0, 10, 0);
  bs._children.push_back(pair<hal_size_t, bool>(0, true));
  bs.applyTo(bi);
     
  ti = child1->getTopSegmentIterator();
  ts.set(0, 10, 0, true, 0);
  ts.applyTo(ti);

  bi = child1->getBottomSegmentIterator();
  bs.set(0, 5, 0);
  bs._children.clear();
  bs.applyTo(bi);
  bi->toRight();
  bs.set(5, 5, 0);
  bs.applyTo(bi);
}

void TopSegmentIteratorReverseTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti, ti2;

  const Genome* parent1 = alignment->openGenome("parent1");
  const Genome* child1 = alignment->openGenome("child1");

  ti = child1->getTopSegmentIterator();
  bi = parent1->getBottomSegmentIterator();

  ti2 = child1->getTopSegmentIterator();
  ti2->toChild(bi, 0);
  
  CuAssertTrue(_testCase, ti->getStartPosition() == 0);
  CuAssertTrue(_testCase, ti->getLength() == 10);
  CuAssertTrue(_testCase, ti->getReversed() == false);

  CuAssertTrue(_testCase, ti2->getStartPosition() == 9);
  CuAssertTrue(_testCase, ti2->getLength() == 10);
  CuAssertTrue(_testCase, ti2->getReversed() == true);

  bi->slice(1, 3);
  ti2->toChild(bi, 0);
  
  CuAssertTrue(_testCase, bi->getStartPosition() == 1);
  CuAssertTrue(_testCase, bi->getLength() == 6);
  CuAssertTrue(_testCase, ti2->getStartPosition() == 8);
  CuAssertTrue(_testCase, ti2->getLength() == 6);

  string buffer;
  bi->getString(buffer);
  CuAssertTrue(_testCase, buffer == "CCTACG");
  ti2->getString(buffer);
  CuAssertTrue(_testCase, buffer == "CACGTA");

  bi = child1->getBottomSegmentIterator();
  CuAssertTrue(_testCase, bi->getReversed() == false);

  ti->toParseUp(bi);  
  CuAssertTrue(_testCase, ti->getStartPosition() == 0);
  CuAssertTrue(_testCase, ti->getLength() == 5);

  bi->toReverse();
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, ti->getStartPosition() == 4);
  CuAssertTrue(_testCase, ti->getLength() == 5);

  bi->toReverse();  
  CuAssertTrue(_testCase, bi->getReversed() == false);
  bi->toRight();
  ti->toParseUp(bi);
  CuAssertTrue(_testCase, ti->getStartPosition() == 5);
  CuAssertTrue(_testCase, ti->getLength() == 5);

  bi->toReverse();  
  ti->toParseUp(bi);  
  CuAssertTrue(_testCase, ti->getStartPosition() == 9);
  CuAssertTrue(_testCase, ti->getLength() == 5);
}

void TopSegmentIsGapTest::createCallBack(AlignmentPtr alignment)
{
  size_t numSequences = 3;
  vector<Sequence::Info> seqVec(numSequences);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  Genome* parent1 = alignment->addRootGenome("parent1");
  Genome* child1 = alignment->addLeafGenome("child1", "parent1", 1);

  // set up two genomes.  each with three sequences.  each sequence
  // with 5 segments of length two.  start with segment i in parent
  // aligned with segment i in child.
  for (size_t i = 0; i < numSequences; ++i)
  {
    stringstream ss;
    ss << "Sequence" << i;
    string name = ss.str();
    seqVec[i] = Sequence::Info(name, 10, 5, 5);
  }
  parent1->setDimensions(seqVec);
  child1->setDimensions(seqVec);

  bi = parent1->getBottomSegmentIterator();
  for (; bi != parent1->getBottomSegmentEndIterator(); bi->toRight())
  {
    bs.set(bi->getBottomSegment()->getArrayIndex() * 2, 2);
    bs._children.clear();
    bs._children.push_back(pair<hal_size_t, bool>(
                            bi->getBottomSegment()->getArrayIndex(), 
                            false));
    bs.applyTo(bi);
  }
     
  ti = child1->getTopSegmentIterator();
  for (; ti != child1->getTopSegmentEndIterator(); ti->toRight())
  {
    ts.set(ti->getTopSegment()->getArrayIndex() * 2, 2, 
           ti->getTopSegment()->getArrayIndex());
    ts.applyTo(ti);
  }

  // insertion in middle (8th top segment)

  bi = parent1->getBottomSegmentIterator(8);
  ti = child1->getTopSegmentIterator(8);
  assert(bi->getBottomSegment()->getChildIndex(0) == 8 &&
         ti->getTopSegment()->getParentIndex() == 8);
  bi->getBottomSegment()->setChildIndex(0, 9);
  ti->getTopSegment()->setParentIndex(NULL_INDEX);
  ti->toRight();
  ti->getTopSegment()->setParentIndex(8);
  
  // insertion at begining (10th top segment)
  
  bi = parent1->getBottomSegmentIterator(10);
  ti = child1->getTopSegmentIterator(10);
  assert(bi->getBottomSegment()->getChildIndex(0) == 10 &&
         ti->getTopSegment()->getParentIndex() == 10);
  bi->getBottomSegment()->setChildIndex(0, 11);
  ti->getTopSegment()->setParentIndex(NULL_INDEX);
  ti->toRight();
  ti->getTopSegment()->setParentIndex(10);

  // just having a null parent is not enough for an insertion
  bi = parent1->getBottomSegmentIterator(2);
  ti = child1->getTopSegmentIterator(2);
  assert(bi->getBottomSegment()->getChildIndex(0) == 2 &&
         ti->getTopSegment()->getParentIndex() == 2);
  ti->getTopSegment()->setParentIndex(NULL_INDEX);
}

void TopSegmentIsGapTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  const Genome* child1 = alignment->openGenome("child1");

  for (hal_size_t i = 0; i < child1->getNumTopSegments(); ++i)
  {
    ti = child1->getTopSegmentIterator(i);
    if (i == 8 || i == 10)
    {
//      CuAssertTrue(_testCase, ti->getTopSegment()->isGapInsertion());
    }
    else
    {
//      CuAssertTrue(_testCase, !ti->getTopSegment()->isGapInsertion());
    }
  }
}

void halTopSegmentSimpleIteratorTest(CuTest *testCase)
{
  try 
  {
    TopSegmentSimpleIteratorTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halTopSegmentSequenceTest(CuTest *testCase)
{
  try 
  {
    TopSegmentSequenceTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halTopSegmentIteratorParseTest(CuTest *testCase)
{
  try 
  {
    TopSegmentIteratorParseTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halTopSegmentIteratorToSiteTest(CuTest *testCase)
{
  try 
  {
    TopSegmentIteratorToSiteTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halTopSegmentIteratorReverseTest(CuTest *testCase)
{
  try 
  {
    TopSegmentIteratorReverseTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halTopSegmentIsGapTest(CuTest *testCase)
{
  try 
  {
    TopSegmentIsGapTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halTopSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halTopSegmentSimpleIteratorTest);
  SUITE_ADD_TEST(suite, halTopSegmentSequenceTest);
  SUITE_ADD_TEST(suite, halTopSegmentIteratorParseTest);
  SUITE_ADD_TEST(suite, halTopSegmentIteratorToSiteTest);
  SUITE_ADD_TEST(suite, halTopSegmentIteratorReverseTest);
  SUITE_ADD_TEST(suite, halTopSegmentIsGapTest);
  return suite;
}

