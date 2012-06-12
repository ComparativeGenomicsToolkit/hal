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
#include "halBottomSegmentTest.h"
#include "halTopSegmentTest.h"
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
      ColumnIterator::DNASet::const_iterator dnaIt = i->second.begin();
      for (size_t j = 0; j < i->second.size(); ++j)
      {
        CuAssertTrue(_testCase, i->second.size() == 1);
        CuAssertTrue(_testCase, (*dnaIt)->getArrayIndex() == 
                     (hal_index_t)columnNumber);
        dnaIt++;
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
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = son2->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = dad->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i, i, 0);
    ts.applyTo(ti);

    bi = dad->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildIndex(1, i);
    bs.set(i * segLength, segLength, i, 0);
    bs.applyTo(bi);

    bi = grandpa->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bs.set(i * segLength, segLength);
    bs.applyTo(bi);
  }  
}

void ColumnIteratorDepthTest::checkGenome(const Genome* genome)
{
  assert(genome != NULL);
  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    CuAssertTrue(_testCase, colMap->size() == 4);
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      CuAssertTrue(_testCase, i->second.size() == 1);
      DNAIteratorConstPtr dnaIt = *i->second.begin();
      
      /*  cout << "column=" << columnNumber 
           << " genome=" << dnaIt->getGenome()->getName()
           << " index=" << dnaIt->getArrayIndex() << endl;*/

      CuAssertTrue(_testCase, dnaIt->getArrayIndex() == 
                   (hal_index_t)columnNumber);
    }
    colIterator->toRight();
  }
}

void ColumnIteratorDepthTest::checkCallBack(AlignmentConstPtr alignment)
{
  // validateAlignment(alignment);
  const Genome* genome = alignment->openGenome("grandpa");
  checkGenome(genome);
  genome = alignment->openGenome("dad");
  checkGenome(genome);
  genome = alignment->openGenome("son1");
  checkGenome(genome);
  genome = alignment->openGenome("son2");
  checkGenome(genome);
}

void ColumnIteratorDupTest::createCallBack(AlignmentPtr alignment)
{
  double branchLength = 1e-10;

  alignment->addRootGenome("dad");
  alignment->addLeafGenome("son1", "dad", branchLength);
  alignment->addLeafGenome("son2", "dad", branchLength);
  
  vector<Sequence::Info> dims(1);
  hal_size_t numSegments = 10;
  hal_size_t segLength = 10;
  hal_size_t seqLength = numSegments * segLength;

  Genome* son1 = alignment->openGenome("son1");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son1->setDimensions(dims);
  son1->setString(string(seqLength, 't'));
  
  Genome* son2 = alignment->openGenome("son2");
  dims[0] = Sequence::Info("seq", seqLength, numSegments, 0);
  son2->setDimensions(dims);
  son2->setString(string(seqLength, 'c'));

  Genome* dad = alignment->openGenome("dad");
  dims[0] = Sequence::Info("seq", seqLength, 0, numSegments);
  dad->setDimensions(dims);
  dad->setString(string(seqLength, 'G'));
  
  BottomSegmentIteratorPtr bi;
  BottomSegmentStruct bs;
  TopSegmentIteratorPtr ti;
  TopSegmentStruct ts;

  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    ti = son2->getTopSegmentIterator(i);
    ts.set(i * segLength, segLength, i);
    ts.applyTo(ti);

    bi = dad->getBottomSegmentIterator(i);
    bi->getBottomSegment()->setChildIndex(0, i);
    bi->getBottomSegment()->setChildIndex(1, i);
    bs.set(i * segLength, segLength);
    bs.applyTo(bi);
  }  
  
  // son1 is just one big duplication stemming 
  // from 0th segment in dad
  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    ti = son1->getTopSegmentIterator(i);
    ti->getTopSegment()->setParentIndex(0);
    if (i < numSegments - 1)
    {
      ti->getTopSegment()->setNextParalogyIndex(i+1);
    }
    else
    {
      ti->getTopSegment()->setNextParalogyIndex(0);
    }
    if (i > 0)
    {
      bi = dad->getBottomSegmentIterator(i);
      bi->getBottomSegment()->setChildIndex(0, NULL_INDEX);
    }
  }
  
  // son2 has a duplication between 4 and 8,
  // which both derive from 4 in dad
  ti = son2->getTopSegmentIterator(4);
  ti->getTopSegment()->setNextParalogyIndex(8);
  ti = son2->getTopSegmentIterator(8);
  ti->getTopSegment()->setNextParalogyIndex(4);
  ti->getTopSegment()->setParentIndex(4);
  bi = dad->getBottomSegmentIterator(8);
  bi->getBottomSegment()->setChildIndex(1, NULL_INDEX);
}

void ColumnIteratorDupTest::checkGenome(const Genome* genome)
{
  assert(genome != NULL);
  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  ColumnIteratorConstPtr endIterator = genome->getColumnEndIterator();
  size_t colNumber = 0;
  for (; colIterator < endIterator; colIterator->toRight())
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    // check that all three genomes are in the map
    CuAssertTrue(_testCase, colMap->size() == 3);

    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      DNAIteratorConstPtr dnaIt = *i->second.begin();
      // the first segment (of any genome) should be aligned to 
      // every segment in son1
      if (i->first->getName() == "son1" && genome->getName() != "son1")
      {
        if (colNumber < i->first->getTopSegmentIterator()->getLength())
        {
          CuAssertTrue(_testCase, i->second.size() ==
                       i->first->getTopSegmentIterator()->getLength());
        }
        else
        {
          CuAssertTrue(_testCase, i->second.size() == 0);
        }
      }
      // check the paralogy on son2 
      else if (i->first->getName() == "son2" && genome->getName() == "dad")
      {
        if (colNumber >= 40 && colNumber < 50)
        {
          CuAssertTrue(_testCase, i->second.size() == 2);
        }
        else if (colNumber >= 80 && colNumber < 90)
        {
          CuAssertTrue(_testCase, i->second.size() == 0);
        }
        else
        {
          CuAssertTrue(_testCase, i->second.size() == 1);
        }
      }
    }
    ++colNumber;
  }
}

void ColumnIteratorDupTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome("dad");
  checkGenome(genome);
  genome = alignment->openGenome("son1");
  checkGenome(genome);
  genome = alignment->openGenome("son2");
  checkGenome(genome);
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
  try 
  {
    ColumnIteratorDepthTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  } 
}

void halColumnIteratorDupTest(CuTest *testCase)
{
  // try 
  {
    ColumnIteratorDupTest tester;
    tester.check(testCase);
  }
  // catch (...) 
  {
//    CuAssertTrue(testCase, false);
  } 
}

CuSuite* halColumnIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halColumnIteratorBaseTest);
  SUITE_ADD_TEST(suite, halColumnIteratorDepthTest);
  SUITE_ADD_TEST(suite, halColumnIteratorDupTest);
  return suite;
}

