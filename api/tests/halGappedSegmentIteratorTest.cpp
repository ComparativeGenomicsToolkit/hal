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

void GappedSegmentSimpleIteratorTest::createCallBack(Alignment* alignment)
{
  addIdenticalParentChild(alignment, 2, 100, 5);
  Genome* parent = alignment->openGenome(alignment->getRootName());
  Genome* child = parent->getChild(0);
  TopSegmentIteratorPtr ti = child->getTopSegmentIterator();
  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
  int i = 0;
  while (not ti->atEnd())
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
GappedSegmentSimpleIteratorTest::checkCallBack(const Alignment* alignment)
{
  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");

  GappedTopSegmentIteratorPtr gtsIt = 
     child->getGappedTopSegmentIterator(0, 9999999);
  GappedBottomSegmentIteratorPtr gbsIt = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
  GappedTopSegmentIteratorPtr gtsItRev = 
     child->getGappedTopSegmentIterator(0, 9999999);
  gtsItRev->toReverse();
  GappedBottomSegmentIteratorPtr gbsItRev = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
  gbsItRev->toReverse();

  for (size_t i = 0; i < child->getNumTopSegments(); ++i)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->equals(gtsIt->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    gtsIt->toRight();

    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->equals(gbsIt->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    gbsIt->toRight();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, tsItRev->equals(gtsItRev->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)tsItRev->getTopSegment()->getArrayIndex() == i);
    gtsItRev->toLeft();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, bsItRev->equals(gbsItRev->getRight()));
    CuAssertTrue(_testCase, 
                 (size_t)bsItRev->getBottomSegment()->getArrayIndex() == i);
    gbsItRev->toLeft();
  }

  gtsIt = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 1, 9999999);
  gbsIt = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 1, 0, 9999999);
  gtsItRev = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 1, 9999999);
  gtsItRev->toReverse();
  gbsItRev = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 1, 0, 9999999);
  gbsItRev->toReverse();

  for (hal_index_t i = child->getNumTopSegments() - 1; i >= 0; --i)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->equals(gtsIt->getRight()));
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gtsIt->getReversed() == false);
    gtsIt->toLeft();

    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->equals(gbsIt->getRight()));
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gbsIt->getReversed() == false);
    gbsIt->toLeft();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, tsItRev->equals(gtsItRev->getRight()));
    CuAssertTrue(_testCase, tsItRev->getTopSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gtsItRev->getReversed() == true);
    gtsItRev->toRight();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, bsItRev->equals(gbsItRev->getRight()));
    CuAssertTrue(_testCase, bsItRev->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gbsItRev->getReversed() == true);
    gbsItRev->toRight();
  }

}

void GappedSegmentSimpleIteratorTest2::createCallBack(Alignment* alignment)
{
  addIdenticalParentChild(alignment, 2, 100, 5);
  Genome* parent = alignment->openGenome(alignment->getRootName());
  Genome* child = parent->getChild(0);
  TopSegmentIteratorPtr ti = child->getTopSegmentIterator();
  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
  hal_index_t i = 0;
  bool reversed = true;
  while (not ti->atEnd())
  {
    if (i % 5 == 0)
    {
      reversed = !reversed;
      if (reversed && i < (hal_index_t)(parent->getNumBottomSegments() - 1))
      {
        makeInversion(ti, 5);
      }
    }

    ti->toRight();
    bi->toRight();
    ++i;
  }
}

void 
GappedSegmentSimpleIteratorTest2::checkCallBack(const Alignment* alignment)
{
  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");

  GappedTopSegmentIteratorPtr gtsIt = 
     child->getGappedTopSegmentIterator(0, 9999999);
  GappedBottomSegmentIteratorPtr gbsIt = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
  GappedTopSegmentIteratorPtr gtsItRev = 
     child->getGappedTopSegmentIterator(0, 9999999);
  gtsItRev->toReverse();
  GappedBottomSegmentIteratorPtr gbsItRev = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
  gbsItRev->toReverse();

  for (size_t i = 0; i < child->getNumTopSegments(); i += 5)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    tsIt = gtsIt->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i + 4);


    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    bsIt = gbsIt->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i + 4);

    GappedBottomSegmentIteratorPtr gappedParent = gbsIt->clone();
    gappedParent->toParent(gtsIt);
    if (gappedParent->getReversed())
    {
      gappedParent->toReverse();
    }
    CuAssertTrue(_testCase, gappedParent->equals(gbsIt));
    GappedTopSegmentIteratorPtr gappedChild = gtsIt->clone();
    gappedChild->toChild(gbsIt);
    if (gappedChild->getReversed())
    {
      gappedChild->toReverse();
    }
    CuAssertTrue(_testCase, gappedChild->equals(gtsIt));
    
    gtsIt->toRight();
    gbsIt->toRight();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)tsItRev->getTopSegment()->getArrayIndex() == i + 4);
    tsItRev = gtsItRev->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)tsItRev->getTopSegment()->getArrayIndex() == i);
    gtsItRev->toLeft();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)bsItRev->getBottomSegment()->getArrayIndex() == i+4);
    bsItRev = gbsItRev->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)bsItRev->getBottomSegment()->getArrayIndex() == i);
    gbsItRev->toLeft();
  }

  gtsIt = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 5, 9999999);
  gbsIt = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 5, 0, 9999999);
  gtsItRev = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 5, 9999999);
  gtsItRev->toReverse();
  gbsItRev = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 5, 0, 9999999);
  gbsItRev->toReverse();

  for (hal_index_t i = child->getNumTopSegments() - 1; i >= 0; i -= 5)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i - 4);
    tsIt = gtsIt->getRight();
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gtsIt->getReversed() == false);
    gtsIt->toLeft();

    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i-4);
    bsIt = gbsIt->getRight();
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gbsIt->getReversed() == false);
    gbsIt->toLeft();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, tsItRev->getTopSegment()->getArrayIndex() == i);
    tsItRev = gtsItRev->getRight();
    CuAssertTrue(_testCase, tsItRev->getTopSegment()->getArrayIndex() == i-4);
    CuAssertTrue(_testCase, gtsItRev->getReversed() == true);
    gtsItRev->toRight();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, bsItRev->getBottomSegment()->getArrayIndex() == i);
    bsItRev = gbsItRev->getRight();
    CuAssertTrue(_testCase, bsItRev->getBottomSegment()->getArrayIndex()==i-4);
    CuAssertTrue(_testCase, gbsItRev->getReversed() == true);
    gbsItRev->toRight();
  }
}

void GappedSegmentIteratorIndelTest::createCallBack(Alignment* alignment)
{
  addIdenticalParentChild(alignment, 1, 20, 5);
  Genome* parent = alignment->openGenome(alignment->getRootName());
  Genome* child = parent->getChild(0);
  TopSegmentIteratorPtr ti = child->getTopSegmentIterator();
  BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator();
//  int i = 0;
//  bool reversed = true;

  bi = parent->getBottomSegmentIterator(0);
  makeDelGap(bi);
  bi = parent->getBottomSegmentIterator(3);
  makeDelGap(bi);
/*
  ti = child->getTopSegmentIterator(1);
  makeInsGap(ti);
  ti = child->getTopSegmentIterator(21);
  makeInsGap(ti);
  ti = child->getTopSegmentIterator(28);
  makeInsGap(ti);
*/  
/*  for (size_t i = 0; i < 20; ++i)
  {
    cout << i << ": ";
    bi = parent->getBottomSegmentIterator(i);
    ti = child->getTopSegmentIterator(i);
    cout << "ci=" << bi->getBottomSegment()->getChildIndex(0) 
         << " pi=" << ti->getTopSegment()->getParentIndex() << endl;
         }*/
}

void 
GappedSegmentIteratorIndelTest::checkCallBack(const Alignment* alignment)
{
  const Genome* child = alignment->openGenome("child");
  const Genome* parent = alignment->openGenome("parent");

  GappedTopSegmentIteratorPtr gtsIt = 
     child->getGappedTopSegmentIterator(0, 9999999);

  GappedBottomSegmentIteratorPtr gbsIt = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
  GappedTopSegmentIteratorPtr gtsItRev = 
     child->getGappedTopSegmentIterator(0, 9999999);
  gtsItRev->toReverse();
  GappedBottomSegmentIteratorPtr gbsItRev = 
     parent->getGappedBottomSegmentIterator(0, 0, 9999999);
     gbsItRev->toReverse();

  for (size_t i = 0; i < child->getNumTopSegments(); i += 20)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i);
    tsIt = gtsIt->getRight();

    CuAssertTrue(_testCase, 
                 (size_t)tsIt->getTopSegment()->getArrayIndex() == i + 19);

    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i);
    bsIt = gbsIt->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)bsIt->getBottomSegment()->getArrayIndex() == i + 19);

    GappedBottomSegmentIteratorPtr gappedParent = gbsIt->clone();
    gappedParent->toParent(gtsIt);
    if (gappedParent->getReversed())
    {
      gappedParent->toReverse();
    }
    CuAssertTrue(_testCase,
                 gappedParent->equals(gbsIt));
    GappedTopSegmentIteratorPtr gappedChild = gtsIt->clone();
    gappedChild->toChild(gbsIt);
    if (gappedChild->getReversed())
    {
      gappedChild->toReverse();
    }
    CuAssertTrue(_testCase, gappedChild->equals(gtsIt));
    
    gtsIt->toRight();
    gbsIt->toRight();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)tsItRev->getTopSegment()->getArrayIndex() == i + 19);
    tsItRev = gtsItRev->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)tsItRev->getTopSegment()->getArrayIndex() == i);
    gtsItRev->toLeft();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, 
                 (size_t)bsItRev->getBottomSegment()->getArrayIndex() == i+19);
    bsItRev = gbsItRev->getRight();
    CuAssertTrue(_testCase, 
                 (size_t)bsItRev->getBottomSegment()->getArrayIndex() == i);
    gbsItRev->toLeft();

    }

  gtsIt = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 20, 9999999);
  gbsIt = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 20, 0, 9999999); 
  gtsItRev = child->getGappedTopSegmentIterator(
    child->getNumTopSegments() - 20, 9999999);
  gtsItRev->toReverse();
  gbsItRev = parent->getGappedBottomSegmentIterator(
    child->getNumTopSegments() - 20, 0, 9999999);
  gbsItRev->toReverse();

  for (hal_index_t i = child->getNumTopSegments() - 1; i >= 0; i -= 20)
  {
    TopSegmentIteratorPtr tsIt = gtsIt->getLeft();
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i - 19);
    tsIt = gtsIt->getRight();
    CuAssertTrue(_testCase, tsIt->getTopSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gtsIt->getReversed() == false);
    gtsIt->toLeft();

    BottomSegmentIteratorPtr bsIt = gbsIt->getLeft();
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i-19);
    bsIt = gbsIt->getRight();
    CuAssertTrue(_testCase, bsIt->getBottomSegment()->getArrayIndex() == i);
    CuAssertTrue(_testCase, gbsIt->getReversed() == false);
    gbsIt->toLeft();

    TopSegmentIteratorPtr tsItRev = gtsItRev->getLeft();
    CuAssertTrue(_testCase, tsItRev->getTopSegment()->getArrayIndex() == i);
    tsItRev = gtsItRev->getRight();
    CuAssertTrue(_testCase, tsItRev->getTopSegment()->getArrayIndex() == i-19);
    CuAssertTrue(_testCase, gtsItRev->getReversed() == true);
    gtsItRev->toRight();

    BottomSegmentIteratorPtr bsItRev = gbsItRev->getLeft();
    CuAssertTrue(_testCase, bsItRev->getBottomSegment()->getArrayIndex() == i);
    bsItRev = gbsItRev->getRight();
    CuAssertTrue(_testCase, bsItRev->getBottomSegment()->getArrayIndex()==i-19);
    CuAssertTrue(_testCase, gbsItRev->getReversed() == true);
    gbsItRev->toRight();
    }

}

void halGappedSegmentSimpleIteratorTest(CuTest *testCase)
{
    GappedSegmentSimpleIteratorTest tester;
    tester.check(testCase);
}

void halGappedSegmentSimpleIteratorTest2(CuTest *testCase)
{
    GappedSegmentSimpleIteratorTest2 tester;
    tester.check(testCase);
}

void halGappedSegmentIteratorIndelTest(CuTest *testCase)
{
    GappedSegmentIteratorIndelTest tester;
    tester.check(testCase);
}

CuSuite* halGappedSegmentIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halGappedSegmentSimpleIteratorTest);
  SUITE_ADD_TEST(suite, halGappedSegmentSimpleIteratorTest2);
  SUITE_ADD_TEST(suite, halGappedSegmentIteratorIndelTest);
  return suite;
}

