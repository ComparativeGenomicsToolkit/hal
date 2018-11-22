/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <iostream>
#include "halAlignmentTest.h"
#include "halValidateTest.h"
#include "halRandomData.h"
#include "halGappedTopSegmentIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halRandNumberGen.h"

extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

static RandNumberGen rng;


void ValidateSmallTest::createCallBack(Alignment* alignment)
{
    createRandomAlignment(rng, alignment, 0.75, 0.1, 2, 5, 10, 1000, 5, 10);
                        
}

void ValidateSmallTest::checkCallBack(const Alignment* alignment)
{
  validateAlignment(alignment);
}

void ValidateMediumTest::createCallBack(Alignment* alignment)
{
    createRandomAlignment(rng, alignment, 1.25, 0.7, 10, 20, 2, 50, 1000, 50000);
                        
}

void ValidateMediumTest::checkCallBack(const Alignment* alignment)
{
  validateAlignment(alignment);
}

void ValidateLargeTest::createCallBack(Alignment* alignment)
{
    createRandomAlignment(rng, alignment, 2.0, 1.0, 50, 100, 2, 10, 10000,
                          500000);
                        
}

void ValidateLargeTest::checkCallBack(const Alignment* alignment)
{
  validateAlignment(alignment);
}

void halValidateSmallTest(CuTest *testCase)
{
    ValidateSmallTest tester;
    tester.check(testCase);
}

void halValidateMediumTest(CuTest *testCase)
{
    ValidateMediumTest tester;
    tester.check(testCase);
}

void halValidateLargeTest(CuTest *testCase)
{
    ValidateLargeTest tester;
    tester.check(testCase);
}

CuSuite* halValidateTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halValidateSmallTest);
  SUITE_ADD_TEST(suite, halValidateMediumTest);
#if 0  // this is very slow
  SUITE_ADD_TEST(suite, halValidateLargeTest);
#endif
  return suite;
}

