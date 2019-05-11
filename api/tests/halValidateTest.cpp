/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halValidateTest.h"
#include "halAlignmentTest.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halRandNumberGen.h"
#include "halRandomData.h"
#include <iostream>
#include <string>

extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

static RandNumberGen rng;

void ValidateSmallTest::createCallBack(Alignment *alignment) {
    createRandomAlignment(rng, alignment, 0.75, 0.1, 2, 5, 10, 1000, 5, 10);
}

void ValidateSmallTest::checkCallBack(const Alignment *alignment) {
    validateAlignment(alignment);
}

void ValidateMediumTest::createCallBack(Alignment *alignment) {
    createRandomAlignment(rng, alignment, 1.25, 0.7, 10, 20, 2, 50, 1000, 50000);
}

void ValidateMediumTest::checkCallBack(const Alignment *alignment) {
    validateAlignment(alignment);
}

void ValidateLargeTest::createCallBack(Alignment *alignment) {
    createRandomAlignment(rng, alignment, 2.0, 1.0, 50, 100, 2, 10, 10000, 500000);
}

void ValidateLargeTest::checkCallBack(const Alignment *alignment) {
    validateAlignment(alignment);
}

void ValidateManyGenomesTest::createCallBack(Alignment *alignment) {
    // hdf5 runs out of memory on Travis.
    if (alignment->getStorageFormat() == "mmap") {
        cout << "running" << endl;
        createRandomAlignment(rng, alignment, 0.75, 0.1, 363, 400, 1, 10, 1, 2);
    }
}

void ValidateManyGenomesTest::checkCallBack(const Alignment *alignment) {
    if (alignment->getStorageFormat() == "mmap") {
        validateAlignment(alignment);
    }
}

void halValidateSmallTest(CuTest *testCase) {
    ValidateSmallTest tester;
    tester.check(testCase);
}

void halValidateMediumTest(CuTest *testCase) {
    ValidateMediumTest tester;
    tester.check(testCase);
}

void halValidateLargeTest(CuTest *testCase) {
    ValidateLargeTest tester;
    tester.check(testCase);
}

void halValidateManyGenomesTest(CuTest *testCase) {
    ValidateManyGenomesTest tester;
    tester.check(testCase);
}

CuSuite *halValidateTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halValidateSmallTest);
    SUITE_ADD_TEST(suite, halValidateMediumTest);
    SUITE_ADD_TEST(suite, halValidateManyGenomesTest);
#if 0 // this is very slow
  SUITE_ADD_TEST(suite, halValidateLargeTest);
#endif
    return suite;
}
