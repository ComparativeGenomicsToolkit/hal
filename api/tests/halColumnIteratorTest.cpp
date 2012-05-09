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

using namespace std;
using namespace hal;

void ColumnIteratorBaseTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 2, 1, 2, 10, 10, 2, 2);
}

void ColumnIteratorBaseTest::checkCallBack(AlignmentConstPtr alignment)
{
  const Genome* genome = alignment->openGenome(alignment->getRootName());
  assert(genome != NULL);

  cout << "genome len " << genome->getSequenceLength() 
       << " num top " << genome->getNumTopSegments()
       << " num bottom " << genome->getNumBottomSegments() << endl;

  ColumnIteratorConstPtr colIterator = genome->getColumnIterator();
  const ColumnIterator::ColumnMap* colMap = colIterator->getColumnMap();
  CuAssertTrue(_testCase, colMap->size() == 2);
  for (ColumnIterator::ColumnMap::const_iterator i = colMap->begin();
       i != colMap->end(); ++i)
  {
    cout << i->second.size() << endl;
    CuAssertTrue(_testCase, i->second.size() == 1);
    CuAssertTrue(_testCase, i->second[0]->getArrayIndex() == 0);
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

