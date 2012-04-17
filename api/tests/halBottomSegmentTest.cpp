/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "halBottomSegmentTest.h"

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
  for (hal_size_t i = i = 0; i < _children.size(); ++i)
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
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    CuAssertTrue(_testCase, tsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, tsIt->getBottomSegment()->getNumChildren() == 
                 numChildren);
    CuAssertTrue(_testCase, _bottomSegments[i]._children.size() == numChildren);
    _bottomSegments[i].applyTo(tsIt);
    tsIt->toRight();
  }
}

void BottomSegmentSimpleIteratorTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openConstGenome("Anc0");
  CuAssertTrue(_testCase, 
               ancGenome->getNumBottomSegments() == _bottomSegments.size());
  BottomSegmentIteratorConstPtr tsIt = ancGenome->getBottomSegmentIterator(0);
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    CuAssertTrue(_testCase, tsIt->getBottomSegment()->getArrayIndex() == i);
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
  const Genome* ancGenome = alignment->openConstGenome("Anc0");
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

CuSuite* halBottomSegmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halBottomSegmentSimpleIteratorTest);
  SUITE_ADD_TEST(suite, halBottomSegmentSequenceTest);
  return suite;
}

