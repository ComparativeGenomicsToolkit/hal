/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include <sstream>
#include "halAlignmentTest.h"
#include "halValidateTest.h"
#include "halRandomData.h"

extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;


void ValidateSmallTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        0.75, 
                        0.1,
                        5,
                        10,
                        1000,
                        5,
                        10);
                        
}

void ValidateSmallTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
}

void ValidateMediumTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        1.25, 
                        0.7,
                        20,
                        2,
                        50,
                        1000,
                        50000);
                        
}

void ValidateMediumTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
}

void ValidateLargeTest::createCallBack(AlignmentPtr alignment)
{
  createRandomAlignment(alignment, 
                        2, 
                        1,
                        100,
                        2,
                        10,
                        10000,
                        500000);
                        
}

void ValidateLargeTest::checkCallBack(AlignmentConstPtr alignment)
{
  validateAlignment(alignment);
}

void halValidateSmallTest(CuTest *testCase)
{
  try
  {
    ValidateSmallTest tester;
    tester.check(testCase);
  }
  catch (hal_exception& e)
  {
    cerr << e.what() << endl;
    CuAssertTrue(testCase, false);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

void halValidateMediumTest(CuTest *testCase)
{
  try
  {
    ValidateMediumTest tester;
    tester.check(testCase);
  }
  catch (hal_exception& e)
  {
    cerr << e.what() << endl;
    CuAssertTrue(testCase, false);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

void halValidateLargeTest(CuTest *testCase)
{
  try
  {
    ValidateLargeTest tester;
    tester.check(testCase);
  }
  catch (hal_exception& e)
  {
    cerr << e.what() << endl;
    CuAssertTrue(testCase, false);
  }
  catch (...) 
  {
    CuAssertTrue(testCase, false);
  }
}

CuSuite* halValidateTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halValidateSmallTest);
  SUITE_ADD_TEST(suite, halValidateMediumTest);
//  SUITE_ADD_TEST(suite, halValidateLargeTest);
  return suite;
}

