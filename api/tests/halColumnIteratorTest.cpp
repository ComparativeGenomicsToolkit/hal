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
  createRandomAlignment(alignment, 2, 1e-10, 2, 10, 10, 3, 3);
}

void ColumnIteratorBaseTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
  const Genome* genome = alignment->openGenome(alignment->getRootName());
  assert(genome != NULL);

  cout << "al size " << alignment->getNumGenomes() << endl;
  cout << alignment->getChildNames(genome->getName())[0] << endl;

  cout << "genome len " << genome->getSequenceLength(); 
  cout << " num top " << genome->getNumTopSegments();
  cout << " num bottom " << genome->getNumBottomSegments();
  cout << " num child " << genome->getNumChildren(); 
  cout << " child len " << genome->getChild(0)->getSequenceLength() << endl;


  cout << "root genome looks like\n";
  BottomSegmentIteratorConstPtr bit = genome->getBottomSegmentIterator();
  cout << "child idx 0= " << bit->getBottomSegment()->getChildIndex(0) << endl;
  bit->toRight();
  cout << "child idx 1= " << bit->getBottomSegment()->getChildIndex(0) << endl;
  bit->toRight();
  cout << "child idx 2= " << bit->getBottomSegment()->getChildIndex(0) << endl;


  cout << "child genome looks like\n";
  const Genome* childGenome = genome->getChild(0);
  TopSegmentIteratorConstPtr tit = childGenome->getTopSegmentIterator();
  cout << "pare idx 0= " << tit->getTopSegment()->getParentIndex() << endl;
  tit->toRight();
  cout << "pare idx 1= " << tit->getTopSegment()->getParentIndex() << endl;
  tit->toRight();
  cout << "pare idx 2= " << tit->getTopSegment()->getParentIndex() << endl;

  // Iterate over the genome, ensuring that base i aligns to single 
  // other base
  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  for (size_t columnNumber = 0; columnNumber < genome->getSequenceLength(); 
       ++columnNumber)
  {
    const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
    CuAssertTrue(_testCase, colMap->size() == 2);
    for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
         i != colMap->end(); ++i)
    {
      for (size_t j = 0; j < i->second.size(); ++j)
      {
        CuAssertTrue(_testCase, i->second.size() == 1);
        CuAssertTrue(_testCase, i->second[j]->getArrayIndex() == columnNumber);
      }
    }

    colIterator->toRight();
  }
}

void halColumnIteratorBaseTestTest(CuTest *testCase)
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

CuSuite* halColumnIteratorTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halColumnIteratorBaseTestTest);
  return suite;
}

