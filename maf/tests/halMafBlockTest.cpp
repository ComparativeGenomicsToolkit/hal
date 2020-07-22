/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMafBlockTest.h"
#include "halMafBlock.h"

using namespace std;
using namespace hal;

void MafBlockCreateTest::createCallBack(AlignmentPtr alignment) {
}

void MafBlockCreateTest::checkCallBack(AlignmentConstPtr alignment) {
}

void halMafBlockCreateTest(CuTest *testCase) {
    try {
#if 0 // FIXME: test callback are empty
    MafBlockCreateTest tester;
    tester.check(testCase);
#else
        std::cerr << "Warning: halMafBlockCreateTest are not implemented" << std::endl;
#endif
    } catch (...) {
        CuAssertTrue(testCase, false);
    }
}

CuSuite *halMafBlockTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halMafBlockCreateTest);
    return suite;
}
