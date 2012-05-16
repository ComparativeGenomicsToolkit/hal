/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include "halColumnIteratorTest.h"
#include "halRandomData.h"
#include "hal.h"


using namespace std;
using namespace hal;

void ColumnIteratorBaseTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 10, 1e-10, 3, 77, 77, 10, 10);
}

void ColumnIteratorBaseTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome(alignment->getRootName());
  assert(genome != NULL);

  // Iterate over the genome, ensuring that base i aligns to single 
  // other base
  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    CuAssertTrue(_testCase, colMap->size() == 3);
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      for (size_t j = 0; j < i->second.size(); ++j)
      {
        CuAssertTrue(_testCase, i->second.size() == 1);
        CuAssertTrue(_testCase, i->second[j]->getArrayIndex() == 
                     (hal_index_t)columnNumber);
      }
    }

    colIterator->toRight();
  }
}

void ColumnIteratorDepthTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("grandpa");
  alignment->addLeafGenome("dad", "grandpa", branchLength);
  alignment->addLeafGenome("son1", "dad", branchLength);
  alignment->addLeafGenome("son2", "dad", branchLength);
  
  vector<Sequence::Info> dims(1);
  hal_size_t numSegments = 10;
  hal_size_t segLength = 10;
  hal_size_t seqLength = numSegments * segLength;

  Genome* son1 = alignment->openGenome("son1");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son1->setDimensions(dims);
  
  Genome* son2 = alignment->openGenome("son2");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son2->setDimensions(dims);

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, numSegments);
  dad->setDimensions(dims);

  Genome* grandpa = alignment->openGenome("grandpa");
  dims[0] = Sequence::Info("seq", seqLength, 0, numSegments);
  grandpa->setDimensions(dims);
  
  BottomSegmentIteratorPtr bit;
  TopSegmentIteratorPtr tit;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    tit = son1->getTopSegmentIterator(i);
    tit->getTopSegment()->setParentIndex(i);
    tit->getTopSegment()->setStartPosition(i * segLength);
    tit->getTopSegment()->setLength(segLength);
    tit->getTopSegment()->setBottomParseIndex(NULL_INDEX);
    tit->getTopSegment()->setBottomParseOffset(0);

    tit = son2->getTopSegmentIterator(i);
    tit->getTopSegment()->setParentIndex(i);
    tit->getTopSegment()->setStartPosition(i * segLength);
    tit->getTopSegment()->setLength(segLength);
    tit->getTopSegment()->setBottomParseIndex(NULL_INDEX);
    tit->getTopSegment()->setBottomParseOffset(0);

    tit = dad->getTopSegmentIterator(i);
    tit->getTopSegment()->setParentIndex(i);
    tit->getTopSegment()->setStartPosition(i * segLength);
    tit->getTopSegment()->setLength(segLength);
    tit->getTopSegment()->setBottomParseIndex(i);
    tit->getTopSegment()->setBottomParseOffset(0);

    bit = dad->getBottomSegmentIterator(i);
    bit->getBottomSegment()->setChildIndex(0, i);
    bit->getBottomSegment()->setChildIndex(1, i);
    bit->getBottomSegment()->setStartPosition(i * segLength);
    bit->getBottomSegment()->setLength(segLength);
    bit->getBottomSegment()->setTopParseIndex(i);
    bit->getBottomSegment()->setTopParseOffset(0);

    bit = grandpa->getBottomSegmentIterator(i);
    bit->getBottomSegment()->setChildIndex(0, i);
    bit->getBottomSegment()->setStartPosition(i * segLength);
    bit->getBottomSegment()->setLength(segLength);
    bit->getBottomSegment()->setTopParseIndex(NULL_INDEX);
    bit->getBottomSegment()->setTopParseOffset(0);
  }  
}

void ColumnIteratorDepthTest::checkCallBack(AlignmentConstPtr alignment)
{
  // validateAlignment(alignment);
  const Genome* genome = alignment->openGenome(alignment->getRootName());
  assert(genome != NULL);

  // Iterate over the root, ensuring that base i aligns to single 
  // other base
  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      CuAssertTrue(_testCase, i->second.size() == 1);
      DNAIteratorConstPtr dnaIt = i->second[0];
      
      cout << "column=" << columnNumber 
           << " genome=" << dnaIt->getGenome()->getName()
           << " index=" << dnaIt->getArrayIndex() << endl;

      CuAssertTrue(_testCase, dnaIt->getArrayIndex() == 
                   (hal_index_t)columnNumber);
    }

    colIterator->toRight();
  }
}

void halColumnIteratorBaseTest(CuTest *testCase)
{
  try 
  {
    ColumnIteratorBaseTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorDepthTest(CuTest *testCase)
{
  // try 
  {
    ColumnIteratorDepthTest tester;
    tester.check(testCase);
  }
  //catch (...) 
  {
    //   CuAssertTrue(testCase, false);
  } 
}

CuSuite* halColumnIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halColumnIteratorBaseTest);
  SUITE_ADD_TEST(suite, halColumnIteratorDepthTest);
  return suite;
}

