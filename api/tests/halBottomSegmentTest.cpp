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
   _arrayIndex = rand();
  _topParseIndex = rand();
  _children.clear();
  for (hal_size_t i = 0; i < numChildren; ++i)
  {
    pair<hal_index_t, bool> child;
    child.first = (bool) rand() % 2;
    child.second = rand();
    _children.push_back(child);
  }
}

void BottomSegmentStruct::set(hal_index_t startPosition,
                              hal_size_t length,
                              hal_index_t topParseIndex)
{
  _startPosition = startPosition;
  _length = length;
  _topParseIndex = topParseIndex;
}

void BottomSegmentStruct::applyTo(BottomSegmentIteratorPtr it) const
{
  BottomSegment* seg = it->getBottomSegment();
  seg->setCoordinates(_startPosition, _length);
  seg->setTopParseIndex(_topParseIndex);
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
  CuAssertTrue(testCase, _topParseIndex == seg->getTopParseIndex());
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
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 10000);
  ancGenome->setDimensions(seqVec);
  
  CuAssertTrue(_testCase, ancGenome->getNumChildren() == numChildren);
  
  _bottomSegments.clear();
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    BottomSegmentStruct bottomSeg;
    bottomSeg.setRandom(numChildren);
    bottomSeg._length = 
       ancGenome->getSequenceLength() / ancGenome->getNumBottomSegments();
    bottomSeg._startPosition = i * bottomSeg._length;
    _bottomSegments.push_back(bottomSeg);
  }
  
  BottomSegmentIteratorPtr bsIt = ancGenome->getBottomSegmentIterator(0);
  BottomSegmentIteratorConstPtr bsEnd = 
     ancGenome->getBottomSegmentEndIterator();
  for (size_t i = 0; bsIt != bsEnd; bsIt->toRight(), ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getNumChildren() == 
                 numChildren);
    CuAssertTrue(_testCase, _bottomSegments[i]._children.size() == numChildren);
    _bottomSegments[i].applyTo(bsIt);
  }
}

void BottomSegmentSimpleIteratorTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  CuAssertTrue(_testCase, 
               ancGenome->getNumBottomSegments() == _bottomSegments.size());
  BottomSegmentIteratorConstPtr bsIt = ancGenome->getBottomSegmentIterator(0);
  for (size_t i = 0; i < ancGenome->getNumBottomSegments(); ++i)
  {
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    _bottomSegments[i].compareTo(bsIt, _testCase);
    bsIt->toRight();
  }
  bsIt = ancGenome->getBottomSegmentIterator(
    ancGenome->getNumBottomSegments() - 1);
  for (hal_index_t i = ancGenome->getNumBottomSegments() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i);
    _bottomSegments[i].compareTo(bsIt, _testCase);
    bsIt->toLeft();
  }

  bsIt = ancGenome->getBottomSegmentIterator(0); 
  bsIt->slice(0, bsIt->getLength() - 1);
  for (hal_index_t i = 0; i < (hal_index_t)ancGenome->getSequenceLength(); ++i)
  {
    CuAssertTrue(_testCase, bsIt->getLength() == 1);
    CuAssertTrue(_testCase, bsIt->getStartPosition() == i);
    bsIt->toRight(bsIt->getStartPosition() + 1);
  }
  bsIt = ancGenome->getBottomSegmentIterator(
    ancGenome->getNumBottomSegments() - 1);
  bsIt->slice(bsIt->getLength() - 1, 0);
  for (hal_index_t i = ancGenome->getSequenceLength() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, bsIt->getLength() == 1);
    CuAssertTrue(_testCase, bsIt->getStartPosition() == i);
    bsIt->toLeft(bsIt->getStartPosition() - 1);
  }

  bsIt = ancGenome->getBottomSegmentIterator(0); 
  bsIt->toReverse();
  CuAssertTrue(_testCase, bsIt->getReversed() == true);
  bsIt->slice(bsIt->getLength() - 1, 0);
  for (hal_index_t i = 0; i < (hal_index_t)ancGenome->getSequenceLength(); ++i)
  {
    CuAssertTrue(_testCase, bsIt->getLength() == 1);
    CuAssertTrue(_testCase, bsIt->getStartPosition() == i);
    bsIt->toLeft(bsIt->getStartPosition() + 1);
  }
  bsIt = ancGenome->getBottomSegmentIterator(
    ancGenome->getNumBottomSegments() - 1);
  bsIt->toReverse();
  bsIt->slice(0, bsIt->getLength() - 1);
  for (hal_index_t i = ancGenome->getSequenceLength() - 1; i >= 0; --i)
  {
    CuAssertTrue(_testCase, bsIt->getLength() == 1);
    CuAssertTrue(_testCase, bsIt->getStartPosition() == i);
    bsIt->toRight(bsIt->getStartPosition() - 1);
  }
}

void BottomSegmentSequenceTest::createCallBack(AlignmentPtr alignment)
{
  Genome* ancGenome = alignment->addRootGenome("Anc0", 0);
  vector<Sequence::Info> seqVec(1);
  seqVec[0] = Sequence::Info("Sequence", 1000000, 5000, 700000);
  ancGenome->setDimensions(seqVec);

  ancGenome->setSubString("CACACATTC", 500, 9);
  BottomSegmentIteratorPtr bsIt = ancGenome->getBottomSegmentIterator(100);
  bsIt->getBottomSegment()->setCoordinates(500, 9);
}

void BottomSegmentSequenceTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* ancGenome = alignment->openGenome("Anc0");
  BottomSegmentIteratorConstPtr bsIt = ancGenome->getBottomSegmentIterator(100);
  CuAssertTrue(_testCase, bsIt->getBottomSegment()->getStartPosition() == 500);
  CuAssertTrue(_testCase, bsIt->getBottomSegment()->getLength() == 9);
  string seq;
  bsIt->getString(seq);
  CuAssertTrue(_testCase, seq == "CACACATTC");
  bsIt->toReverse();
  bsIt->getString(seq);
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
  ts.set(0, 10, NULL_INDEX, false, 0);
  ts.applyTo(ti);
  
  bi = case1->getBottomSegmentIterator();
  bs.set(0, 10, 0);
  bs.applyTo(bi);

  // case 2: bottom segment is completely contained in top segment
  Genome* case2 = alignment->addRootGenome("case2");
  seqVec[0] = Sequence::Info("Sequence", 10, 2, 3);
  case2->setDimensions(seqVec);
  
  ti = case2->getTopSegmentIterator();
  ts.set(0, 9, NULL_INDEX, false, 0);
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

void BottomSegmentIteratorToSiteTest::createCallBack(AlignmentPtr alignment)
{
  vector<Sequence::Info> seqVec(1);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  
  // case 1: single segment
  Genome* case1 = alignment->addRootGenome("case1");
  seqVec[0] = Sequence::Info("Sequence", 10, 0, 2);
  case1->setDimensions(seqVec);
  bi = case1->getBottomSegmentIterator();
  bs.set(0, 9);
  bs.applyTo(bi);
  bi->toRight();
  bs.set(9, 1);
  bs.applyTo(bi);
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
  seqVec[0] = Sequence::Info("Sequence", total, 0, numSegs);
  case2->setDimensions(seqVec);
  hal_index_t prev = 0;
  for (size_t i = 0 ; i < numSegs; ++i)
  {
    bi = case2->getBottomSegmentIterator((hal_index_t)i);
    bs.set(prev, segLens[i]);
    prev += segLens[i];
    bs.applyTo(bi);
  }
}

void BottomSegmentIteratorToSiteTest::checkGenome(const Genome* genome)
{
  BottomSegmentIteratorConstPtr bi = genome->getBottomSegmentIterator();
  for (hal_index_t pos = 0; 
       pos < (hal_index_t)genome->getSequenceLength(); ++pos)
  {
    bi->toSite(pos);
    CuAssertTrue(_testCase, bi->getStartPosition() == pos);
    CuAssertTrue(_testCase, bi->getLength() == 1);
    bi->toSite(pos, false);
    CuAssertTrue(_testCase, pos >= bi->getStartPosition() && 
                 pos < bi->getStartPosition() + (hal_index_t)bi->getLength());
    CuAssertTrue(_testCase, 
                 bi->getLength() == bi->getBottomSegment()->getLength());
  }
}

void BottomSegmentIteratorToSiteTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;

  // case 1
  const Genome* case1 = alignment->openGenome("case1");
  checkGenome(case1);

  // case 2
  const Genome* case2 = alignment->openGenome("case2");
  checkGenome(case2);


}

void BottomSegmentIteratorReverseTest::createCallBack(AlignmentPtr alignment)
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

  ti = parent1->getTopSegmentIterator();
  ts.set(0, 5, 0, false, 0);
  ts.applyTo(ti);
  ti->toRight();
  ts.set(5, 5, 0, false, 0);
  ts.applyTo(ti);
}

void BottomSegmentIteratorReverseTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi, bi2;
  TopSegmentIteratorConstPtr ti;

  const Genome* parent1 = alignment->openGenome("parent1");
  const Genome* child1 = alignment->openGenome("child1");

  ti = child1->getTopSegmentIterator();
  bi = parent1->getBottomSegmentIterator();

  bi2 = parent1->getBottomSegmentIterator();
  bi2->toParent(ti);
  
  CuAssertTrue(_testCase, bi->getStartPosition() == 0);
  CuAssertTrue(_testCase, bi->getLength() == 10);
  CuAssertTrue(_testCase, bi->getReversed() == false);

  CuAssertTrue(_testCase, bi2->getStartPosition() == 9);
  CuAssertTrue(_testCase, bi2->getLength() == 10);
  CuAssertTrue(_testCase, bi2->getReversed() == true);

  ti->slice(1, 3);
  bi2->toParent(ti);
  
  CuAssertTrue(_testCase, ti->getStartPosition() == 1);
  CuAssertTrue(_testCase, ti->getLength() == 6);
  CuAssertTrue(_testCase, bi2->getStartPosition() == 8);
  CuAssertTrue(_testCase, bi2->getLength() == 6);

  string buffer;
  ti->getString(buffer);
  CuAssertTrue(_testCase, buffer == "CCTACG");
  bi2->getString(buffer);
  CuAssertTrue(_testCase, buffer == "CACGTA");

  ti = parent1->getTopSegmentIterator();
  CuAssertTrue(_testCase, ti->getReversed() == false);

  bi->toParseDown(ti);  
  CuAssertTrue(_testCase, bi->getStartPosition() == 0);
  CuAssertTrue(_testCase, bi->getLength() == 5);

  ti->toReverse();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == 4);
  CuAssertTrue(_testCase, bi->getLength() == 5);

  ti->toReverse();  
  CuAssertTrue(_testCase, ti->getReversed() == false);
  ti->toRight();
  bi->toParseDown(ti);
  CuAssertTrue(_testCase, bi->getStartPosition() == 5);
  CuAssertTrue(_testCase, bi->getLength() == 5);

  ti->toReverse();  
  bi->toParseDown(ti);  
  CuAssertTrue(_testCase, bi->getStartPosition() == 9);
  CuAssertTrue(_testCase, bi->getLength() == 5);
}

void BottomSegmentIsGapTest::createCallBack(AlignmentPtr alignment)
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

  // deletion in middle (3rd bottom segment)

  bi = parent1->getBottomSegmentIterator(3);
  ti = child1->getTopSegmentIterator(3);
  assert(bi->getBottomSegment()->getChildIndex(0) == 3 &&
         ti->getTopSegment()->getParentIndex() == 3);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  ti->getTopSegment()->setParentIndex(4);
  bi->toRight();
  bi->getBottomSegment()->setChildIndex(0, 3);

  // deletion at beginning (5th bottom segment)
  
  bi = parent1->getBottomSegmentIterator(5);
  ti = child1->getTopSegmentIterator(5);
  assert(bi->getBottomSegment()->getChildIndex(0) == 5 &&
         ti->getTopSegment()->getParentIndex() == 5);
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
  bi->toRight();
  bi->getBottomSegment()->setChildIndex(0, 5);
  ti->toRight();
  ti->getTopSegment()->setParentIndex(5);
  ti->toLeft();
  ti->toLeft();
  ti->getTopSegment()->setParentIndex(NULL_INDEX);
}

void BottomSegmentIsGapTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  const Genome* parent1 = alignment->openGenome("parent1");

  for (hal_size_t i = 0; i < parent1->getNumBottomSegments(); ++i)
  {
    bi = parent1->getBottomSegmentIterator(i);
    if (i == 5 || i == 3)
    {
      // CuAssertTrue(_testCase, bi->getBottomSegment()->isGapDeletion(0));
    }
    else
    {
      // CuAssertTrue(_testCase, !bi->getBottomSegment()->isGapDeletion(0));
    }
  }
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

void halBottomSegmentIteratorToSiteTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentIteratorToSiteTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halBottomSegmentIteratorReverseTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentIteratorReverseTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halBottomSegmentIsGapTest(CuTest *testCase)
{
  try 
  {
    BottomSegmentIsGapTest tester;
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
  SUITE_ADD_TEST(suite, halBottomSegmentIteratorToSiteTest);
  SUITE_ADD_TEST(suite, halBottomSegmentIteratorReverseTest);
  SUITE_ADD_TEST(suite, halBottomSegmentIsGapTest);
  return suite;
}

