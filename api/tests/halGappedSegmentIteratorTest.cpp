/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "halGappedSegmentIteratorTest.h"
#include "halTopSegmentTest.h"
#include "halBottomSegmentTest.h"
#include "halRearrangementTest.h"

using namespace std;
using namespace hal;

void GappedSegmentSimpleIteratorTest::createCallBack(AlignmentPtr alignment)
{
  addIdenticalParentChild(alignment, 2, 100, 5);
  Genome* parent = alignment->openGenome(alignment->getRootName());
  Genome* child = parent->getChild(0);
  TopSegmentIteratorPtr ti = child->getTopSegmentIterator();
  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
  int i = 0;
  while (ti != child->getTopSegmentEndIterator())
  {
    if (i++ % 2)
    {
      ti->getTopSegment()->setParentReversed(true);
      bi->getBottomSegment()->setChildReversed(0, true);
    }
    ti->toRight();
    bi->toRight();
  }
}

void 
GappedSegmentSimpleIteratorTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");
  GappedTopSegmentIteratorConstPtr gtsIt = 
     child->getGappedTopSegmentIterator(0, 100000000000);
  GappedBottomSegmentIteratorConstPtr gbsIt = 
     parent->getGappedBottomSegmentIterator(0, 0, 100000000000);
  for (size_t i = 0; i < child->getNumTopSegments(); ++i)
  {
    TopSegmentIteratorConstPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->equals(gtsIt->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    gtsIt->toRight();

    BottomSegmentIteratorConstPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->equals(gbsIt->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    gbsIt->toRight();
    
  }
  gtsIt = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 1, 100000000000);
  gbsIt = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 1, 0, 100000000000);

  for (hal_index_t i = child->getNumTopSegments() - 1; i >= 0; --i)
  {
    TopSegmentIteratorConstPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->equals(gtsIt->getRight()));
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    gtsIt->toLeft();

    BottomSegmentIteratorConstPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->equals(gbsIt->getRight()));
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i);
    gbsIt->toLeft();

  }

}

void halGappedSegmentSimpleIteratorTest(CuTest *testCase)
{
//  try 
  {
    GappedSegmentSimpleIteratorTest tester;
    tester.check(testCase);
  }
  //catch (...) 
  {
    //CuAssertTrue(testCase, false);
  } 
}

CuSuite* halGappedSegmentIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halGappedSegmentSimpleIteratorTest);
  return suite;
}

