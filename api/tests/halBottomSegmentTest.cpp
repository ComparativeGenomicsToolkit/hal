/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cassert>
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"

using namespace std;
using namespace hal;

// just create a bunch of garbage data.  we don't care 
// about logical consistency for this test, just whether or not
// it's read and written properly. 
void BottomSegmentStruct::setRandom(hal_size_t numChildren)
{
  _length = rand();
  _startPosition = rand();
  _nextParalogyIndex = rand();
   _arrayIndex = rand();
  _topParseIndex = rand();
  _topParseOffset = rand();
  _children.clear();
  for (hal_size_t i = 0; i < numChildren; ++i)
  {
    pair<hal_index_t, hal_bool_t> child;
    child.first = (hal_bool_t) rand() % 2;
    child.second = rand();
    _children.push_back(child);
  }
}

void BottomSegmentStruct::set(hal_index_t startPosition,
                              hal_size_t length,
                              hal_index_t topParseIndex,
                              hal_offset_t topParseOffset,
                              hal_index_t nextParalogyIndex)
{
  _startPosition = startPosition;
  _length = length;
  _topParseIndex = topParseIndex;
  _topParseOffset = topParseOffset;
  _nextParalogyIndex = nextParalogyIndex;
}

void BottomSegmentStruct::applyTo(BottomSegmentIteratorPtr it) const
{
  BottomSegment* seg = it->getBottomSegment();
  seg->setLength(_length);
  seg->setStartPosition(_startPosition);
  seg->setNextParalogyIndex(_nextParalogyIndex);
  seg->setTopParseIndex(_topParseIndex);
  seg->setTopParseOffset(_topParseOffset);
  for (hal_size_t i = 0; i < _children.size(); ++i)
  {
    seg->setChildIndex(i, _children[i].first);
    seg->setChildReversed(i, _children[i].second);
  }
}

void BottomSegmentStruct::compareTo(BottomSegmentIteratorConstPtr it, 
                                 CuTest* testCase) const
{
  const BottomSegment* seg = it->getBottomSegment();
  CuAssertTrue(testCase, _length == seg->getLength());
  CuAssertTrue(testCase, _startPosition == seg->getStartPosition());
  CuAssertTrue(testCase, _nextParalogyIndex == seg->getNextParalogyIndex());
  CuAssertTrue(testCase, _topParseIndex == seg->getTopParseIndex());
  CuAssertTrue(testCase, _topParseOffset == seg->getTopParseOffset());
  CuAssertTrue(testCase, _children.size() == seg->getNumChildren());
  for (hal_size_t i = 0; i < _children.size(); ++i)
  {
    CuAssertTrue(testCase, _children[i].first == seg->getChildIndex(i));
    CuAssertTrue(testCase, _children[i].second == seg->getChildReversed(i));
  }
}

void BottomSegmentSimpleIteratorTest::createCallBack(AlignmentPtr alignment)
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
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);
  
  CuAssertTrue(_testCase, ancGenome->getNumChildren() == numChildren);
  
  _bottomSegments.clear();
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    BottomSegmentStruct bottomSeg;
    bottomSeg.setRandom(numChildren);
    _bottomSegments.push_back(bottomSeg);
  }
  
  BottomSegmentIteratorPtr tsIt = ancGenome->getBottomSegmentIterator(0);
  BottomSegmentIteratorConstPtr tsEnd = 
     ancGenome->getBottomSegmentEndIterator();
  for (size_t i = 0; tsIt != tsEnd; tsIt->toRight(), ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, tsIt->getBottomSegment()->getNumChildren() == 
                 numChildren);
    CuAssertTrue(_testCase, _bottomSegments[i]._children.size() == numChildren);
    _bottomSegments[i].applyTo(tsIt);
  }
}

void BottomSegmentSimpleIteratorTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  CuAssertTrue(_testCase, 
               ancGenome->getNumBottomSegments() == _bottomSegments.size());
  BottomSegmentIteratorConstPtr tsIt = ancGenome->getBottomSegmentIterator(0);
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getBottomSegment()->getArrayIndex() == i);
    _bottomSegments[i].compareTo(tsIt, _testCase);
    tsIt->toRight();
  }
  tsIt = ancGenome->getBottomSegmentIterator(
    ancGenome->getNumBottomSegments() - 1);
  for (hal_index_t i = ancGenome->getNumBottomSegments() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, tsIt->getBottomSegment()->getArrayIndex() == i);
    _bottomSegments[i].compareTo(tsIt, _testCase);
    tsIt->toLeft();
  }
}

void BottomSegmentSequenceTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("Anc0", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);

  ancGenome->setSubString("CACACATTC", 500, 9);
  BottomSegmentIteratorPtr tsIt = ancGenome->getBottomSegmentIterator(100);
  tsIt->getBottomSegment()->setStartPosition(500);
  tsIt->getBottomSegment()->setLength(9);
}

void BottomSegmentSequenceTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  BottomSegmentIteratorConstPtr tsIt = ancGenome->getBottomSegmentIterator(100);
  CuAssertTrue(_testCase, tsIt->getBottomSegment()->getStartPosition() == 500);
  CuAssertTrue(_testCase, tsIt->getBottomSegment()->getLength() == 9);
  string seq;
  tsIt->getString(seq);
  CuAssertTrue(_testCase, seq == "CACACATTC");
  tsIt->toReverse();
  tsIt->getString(seq);
  CuAssertTrue(_testCase, seq == "GAATGTGTG");
}

void BottomSegmentIteratorParseTest::createCallBack(AlignmentPtr alignment)
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
  bs.set(3, 4, 0, 0);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(7, 3, 0, 0);
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
  bs.set(5, 5, 0, 0);
  bs.applyTo(bi);
}

void BottomSegmentIteratorParseTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  // case 1
  const Genome* case1 = alignment->openGenome("case1");
  ti = case1->getTopSegmentIterator();
  bi = case1->getBottomSegmentIterator();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  CuAssertTrue(_testCase, bi->getLength() == ti->getLength());
  ti->slice(3, 1);
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, ti->getLength() == 
               ti->getTopSegment()->getLength() - 4);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  CuAssertTrue(_testCase, bi->getLength() == ti->getLength());

  // case 2
  const Genome* case2 = alignment->openGenome("case2");
  ti = case2->getTopSegmentIterator();
  bi = case2->getBottomSegmentIterator();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  ti->slice(5, 2);
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());

  // case 3
  const Genome* case3 = alignment->openGenome("case3");
  ti = case3->getTopSegmentIterator(1);
  bi = case3->getBottomSegmentIterator();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  ti->slice(2, 1);
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());

  // case 4
 const Genome* case4 = alignment->openGenome("case4");
  ti = case4->getTopSegmentIterator();
  bi = case4->getBottomSegmentIterator();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
  ti->slice(6, 2);
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == ti->getStartPosition());
}

void halBottomSegmentSimpleIteratorTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentSimpleIteratorTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halBottomSegmentSequenceTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentSequenceTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halBottomSegmentIteratorParseTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentIteratorParseTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halBottomSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halBottomSegmentSimpleIteratorTest);
  SUITE_ADD_TEST(suite, halBottomSegmentSequenceTest);
  SUITE_ADD_TEST(suite, halBottomSegmentIteratorParseTest);
  return suite;
}

