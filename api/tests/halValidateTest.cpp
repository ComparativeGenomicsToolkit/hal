/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "halApiTestSupport.h"
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

struct ValidateSmallTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        createRandomAlignment(rng, alignment, 0.75, 0.1, 2, 5, 10, 1000, 5, 10);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        validateAlignment(alignment.get());
    }
};

struct ValidateMediumTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        createRandomAlignment(rng, alignment, 1.25, 0.7, 10, 20, 2, 50, 1000, 50000);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        validateAlignment(alignment.get());
    }
};

struct ValidateLargeTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        createRandomAlignment(rng, alignment, 2.0, 1.0, 50, 100, 2, 10, 10000, 500000);
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        validateAlignment(alignment.get());
    }

};

struct ValidateManyGenomesTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment) {
        // hdf5 runs out of memory on Travis.
        if (alignment->getStorageFormat() == "mmap") {
            cout << "running" << endl;
            createRandomAlignment(rng, alignment, 0.75, 0.1, 363, 400, 1, 10, 1, 2);
        }
    }

    void checkCallBack(AlignmentConstPtr alignment) {
        if (alignment->getStorageFormat() == "mmap") {
            validateAlignment(alignment.get());
        }
    }
};

static void halValidateSmallTest(CuTest *testCase) {
    ValidateSmallTest tester;
    tester.check(testCase);
}

static void halValidateMediumTest(CuTest *testCase) {
    ValidateMediumTest tester;
    tester.check(testCase);
}

static void halValidateLargeTest(CuTest *testCase) {
    ValidateLargeTest tester;
    tester.check(testCase);
}

static void halValidateManyGenomesTest(CuTest *testCase) {
    ValidateManyGenomesTest tester;
    tester.check(testCase);
}

static CuSuite *halValidateTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, halValidateSmallTest);
    SUITE_ADD_TEST(suite, halValidateMediumTest);
    SUITE_ADD_TEST(suite, halValidateManyGenomesTest);
    if (false) {// FIXME: this is very slow
        SUITE_ADD_TEST(suite, halValidateLargeTest);
    }
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halValidateTestSuite());
}
