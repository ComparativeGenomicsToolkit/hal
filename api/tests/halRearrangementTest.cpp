/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "halRearrangementTest.h"
#include "halTopSegmentTest.h"
#include "halBottomSegmentTest.h"

using namespace std;
using namespace hal;

void addIdenticalParentChild(hal::AlignmentPtr alignment,
                             size_t numSequences,
                             size_t numSegmentsPerSequence,
                             size_t segmentLength)
{
  vector<Sequence::Info> seqVec(numSequences);
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;
  
  Genome* parent = alignment->addRootGenome("parent");
  Genome* child = alignment->addLeafGenome("child", "parent", 1);

  for (size_t i = 0; i < numSequences; ++i)
  {
    stringstream ss;
    ss << "Sequence" << i;
    string name = ss.str();
    seqVec[i] = Sequence::Info(name, segmentLength * numSegmentsPerSequence, 
                               numSegmentsPerSequence, 
                               numSegmentsPerSequence);
  }
  parent->setDimensions(seqVec);
  child->setDimensions(seqVec);

  bi = parent->getBottomSegmentIterator();
  for (; bi != parent->getBottomSegmentEndIterator(); bi->toRight())
  {
    bs.set(bi->getBottomSegment()->getArrayIndex() * segmentLength, 
           segmentLength);
    bs._children.clear();
    bs._children.push_back(pair<hal_size_t, bool>(
                            bi->getBottomSegment()->getArrayIndex(), 
                            false));
    bs.applyTo(bi);
  }
     
  ti = child->getTopSegmentIterator();
  for (; ti != child->getTopSegmentEndIterator(); ti->toRight())
  {
    ts.set(ti->getTopSegment()->getArrayIndex() * segmentLength, 
           segmentLength, 
           ti->getTopSegment()->getArrayIndex());
    ts.applyTo(ti);
  } 
}

// doesn't currently work on sequence endpoints
void makeInsertion(BottomSegmentIteratorPtr bi)
{
  assert(bi->getBottomSegment()->isLast() == false);
  TopSegmentIteratorPtr ti =
     bi->getBottomSegment()->getGenome()->getTopSegmentIterator();
  ti->toChild(bi, 0);
  hal_index_t pi = ti->getTopSegment()->getParentIndex();
  ti->getTopSegment()->setParentIndex(NULL_INDEX);
  ti->toRight();
  ti->getTopSegment()->setParentIndex(pi);

  hal_index_t ci = bi->getBottomSegment()->getChildIndex(0);
  bi->getBottomSegment()->setChildIndex(0, ci + 1);
  bi->toRight();
  bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
}

void RearrangementInsertionTest::createCallBack(AlignmentPtr alignment)
{
  size_t numSequences = 3;
  size_t numSegmentsPerSequence = 10;
  size_t segmentLength = 5;
  
  addIdenticalParentChild(alignment, numSequences, numSegmentsPerSequence,
                          segmentLength);

  Genome* parent = alignment->openGenome("parent");
  Genome* child = alignment->openGenome("child");

  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
  
  // insertion smaller than gap threshold
  makeInsertion(bi);

  // insertion larger than gap threshold
  bi->toRight();
  bi->toRight();
  makeInsertion(bi);
  bi->toRight();
  makeInsertion(bi);

  // insertion larger than gap threshold but that contains gaps
}

void RearrangementInsertionTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");
  
  RearrangementPtr r = child->getRearrangement();
//  bool res = r->identifyFromLeftBreakpoint(parent->getTopSegmentIterator());
  
}

void RearrangementSimpleInversionTest::createCallBack(AlignmentPtr alignment)
{
  size_t numSequences = 3;
  size_t numSegmentsPerSequence = 10;
  size_t segmentLength = 50;
  
  addIdenticalParentChild(alignment, numSequences, numSegmentsPerSequence,
                          segmentLength);

  Genome* parent = alignment->openGenome("parent");
  Genome* child = alignment->openGenome("child");

  // inversion on 1st segment (interior)
  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator(1);
  TopSegmentIteratorPtr ti = child->getTopSegmentIterator(1);
  bi->getBottomSegment()->setChildReversed(0, true);
  ti->getTopSegment()->setParentReversed(true);
  
  // inversion on last segment of first sequence
  bi = parent->getBottomSegmentIterator(numSegmentsPerSequence - 1);
  ti = child->getTopSegmentIterator(numSegmentsPerSequence - 1);
  bi->getBottomSegment()->setChildReversed(0, true);
  ti->getTopSegment()->setParentReversed(true);

  // inversion on first segment of 3rd sequence
  bi = parent->getBottomSegmentIterator(numSegmentsPerSequence * 2);
  ti = child->getTopSegmentIterator(numSegmentsPerSequence * 2);
  bi->getBottomSegment()->setChildReversed(0, true);
  ti->getTopSegment()->setParentReversed(true);  
}

void 
RearrangementSimpleInversionTest::checkCallBack(AlignmentConstPtr alignment)
{
  BottomSegmentIteratorConstPtr bi;
  TopSegmentIteratorConstPtr ti;

  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");
  
  size_t numSequences = child->getNumSequences();
  size_t numSegmentsPerSequence = 
     child->getSequenceIterator()->getSequence()->getNumTopSegments();
  size_t segmentLength = 
     child->getTopSegmentIterator()->getTopSegment()->getLength();
  
  RearrangementPtr r = child->getRearrangement();
  do
  {
    hal_index_t leftIdx = 
       r->getLeftBreakpoint()->getTopSegment()->getArrayIndex();
    cout << leftIdx << ") " << r->getID() << endl;
    if (leftIdx == 1 || leftIdx == numSegmentsPerSequence - 1 ||
        leftIdx == (numSegmentsPerSequence * 2) )
    {
      CuAssertTrue(_testCase, r->getID() == Rearrangement::Inversion);
    }
    else
    {
      CuAssertTrue(_testCase, r->getID() != Rearrangement::Inversion);
    }
  }
  while (r->identifyNext() == true);
}

void halRearrangementInsertionTest(CuTest *testCase)
{
  try 
  {
    RearrangementInsertionTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halRearrangementSimpleInversionTest(CuTest *testCase)
{
  try 
  {
    RearrangementSimpleInversionTest tester;
    tester.check(testCase);
  }
   catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halRearrangementTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halRearrangementInsertionTest);
  SUITE_ADD_TEST(suite, halRearrangementSimpleInversionTest);
  return suite;
}

