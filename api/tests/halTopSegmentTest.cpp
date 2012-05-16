/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
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
  _bottomParseOffset = rand();
}

void TopSegmentStruct::set(hal_index_t startPosition,
                           hal_size_t length,
                           hal_index_t parentIndex,
                           hal_index_t bottomParseIndex,
                           hal_offset_t bottomParseOffset,
                           hal_index_t nextParalogyIndex)
{
  _startPosition = startPosition;
  _length = length;
  _parentIndex = parentIndex;
  _bottomParseIndex = bottomParseIndex;
  _bottomParseOffset = bottomParseOffset;
  _nextParalogyIndex = nextParalogyIndex;
}

void TopSegmentStruct::applyTo(TopSegmentIteratorPtr it) const
{
  TopSegment* seg = it->getTopSegment();
  seg->setLength(_length);
  seg->setStartPosition(_startPosition);
  seg->setNextParalogyIndex(_nextParalogyIndex);
  seg->setParentIndex(_parentIndex);
  seg->setBottomParseIndex(_bottomParseIndex);
  seg->setBottomParseOffset(_bottomParseOffset);
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
  CuAssertTrue(testCase, _bottomParseOffset == seg->getBottomParseOffset());
}

void TopSegmentSimpleIteratorTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("Anc0", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);
  
  _topSegments.clear();
  for (size_t i = 0; i < ancGenome->getNumTopSegments(); ++i)
  {
    TopSegmentStruct topSeg;
    topSeg.setRandom();
    _topSegments.push_back(topSeg);
  }
  
  TopSegmentIteratorPtr tsIt = ancGenome->getTopSegmentIterator();
  TopSegmentIteratorConstPtr tsEnd = ancGenome->getTopSegmentEndIterator();
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
  CuAssertTrue(_testCase, ancGenome->getNumTopSegments() == _topSegments.size());
  TopSegmentIteratorConstPtr tsIt = ancGenome->getTopSegmentIterator(0);
  for (size_t i = 0; i < ancGenome->getNumTopSegments(); ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    _topSegments[i].compareTo(tsIt, _testCase);
    tsIt->toRight();
  }
  tsIt = ancGenome->getTopSegmentIterator(ancGenome->getNumTopSegments() - 1);
  for (hal_index_t i = ancGenome->getNumTopSegments() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    _topSegments[i].compareTo(tsIt, _testCase);
    tsIt->toLeft();
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
  tsIt->getTopSegment()->setStartPosition(500);
  tsIt->getTopSegment()->setLength(9);
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
  ts.set(0, 10, NULL_INDEX, 0, 0);
  ts.applyTo(ti);
  
  bi = case1->getBottomSegmentIterator();
  bs.set(0, 10, 0, 0);
  bs.applyTo(bi);

  // case 2: bottom segment is completely contained in top segment
  Genome* case2 = alignment->addRootGenome("case2");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 3);
  case2->setDimensions(seqVec);
  
  ti = case2->getTopSegmentIterator();
  ts.set(0, 9, NULL_INDEX, 0, 0);
  ts.applyTo(ti);

  bi = case2->getBottomSegmentIterator();
  bs.set(0, 3, 0, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(3, 4, 0, 3);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(7, 3, 0, 7);
  bs.applyTo(bi);

  // case 3 top segment is completely contained in bottom segment
  Genome* case3 = alignment->addRootGenome("case3");
  seqVec[0] = Sequence::Info("Sequence", 10, 3, 2);
  case3->setDimensions(seqVec);

  ti = case3->getTopSegmentIterator();
  ts.set(0, 3, NULL_INDEX, 0, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(3, 4, NULL_INDEX, 0, 3);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(7, 3, NULL_INDEX, 0, 7);
  ts.applyTo(ti);

  bi = case3->getBottomSegmentIterator();
  bs.set(0, 9, 0, 0);
  bs.applyTo(bi);
 
  // case 4: top segment overhangs bottom segment on the left
  Genome* case4 = alignment->addRootGenome("case4");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 2);
  case4->setDimensions(seqVec);

  ti = case4->getTopSegmentIterator();
  ts.set(0, 9, NULL_INDEX, 0, 0);
  ts.applyTo(ti);

  bi = case4->getBottomSegmentIterator();
  bs.set(0, 5, 0, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(5, 5, 0, 5);
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

CuSuite* halTopSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halTopSegmentSimpleIteratorTest);
  SUITE_ADD_TEST(suite, halTopSegmentSequenceTest);
  SUITE_ADD_TEST(suite, halTopSegmentIteratorParseTest);
  return suite;
}

