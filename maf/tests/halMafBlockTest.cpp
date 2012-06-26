/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafBlockTest.h"
#include "halMafBlock.h"

using namespace std;
using namespace hal;


void MafBlockCreateTest::createCallBack(AlignmentPtr alignment)
{
  
}

void MafBlockCreateTest::checkCallBack(AlignmentConstPtr alignment)
{
 
}

void halMafBlockCreateTest(CuTest *testCase)
{
  try
  {
    MafBlockCreateTest tester;
    tester.check(testCase);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}


CuSuite* halMafBlockTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halMafBlockCreateTest);
  return suite;
}

