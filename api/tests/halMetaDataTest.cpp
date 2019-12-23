/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halApiTestSupport.h"
#include "halAlignment.h"
#include "halGenome.h"
#include "halMetaData.h"
#include <iostream>
#include <string>
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

struct MetaDataTest : public AlignmentTest {
    void createCallBack(Alignment *alignment) {
        hal_size_t alignmentSize = alignment->getNumGenomes();
        CuAssertTrue(_testCase, alignmentSize == 0);

        MetaData *meta = alignment->getMetaData();
        CuAssertTrue(_testCase, meta->getMap().empty() == true);
        meta->set("colour", "red");
        meta->set("number", "1");
        meta->set("animal", "cat");
        meta->set("colour", "black");

        CuAssertTrue(_testCase, meta->get("colour") == "black");
        CuAssertTrue(_testCase, meta->get("number") == "1");
        CuAssertTrue(_testCase, meta->get("animal") == "cat");

        CuAssertTrue(_testCase, meta->has("colour") == true);
        CuAssertTrue(_testCase, meta->has("city") == false);

        CuAssertTrue(_testCase, meta->getMap().size() == 3);
    }

    void checkCallBack(const Alignment *alignment) {
        const MetaData *meta = alignment->getMetaData();

        CuAssertTrue(_testCase, meta->get("colour") == "black");
        CuAssertTrue(_testCase, meta->get("number") == "1");
        CuAssertTrue(_testCase, meta->get("animal") == "cat");

        CuAssertTrue(_testCase, meta->has("colour") == true);
        CuAssertTrue(_testCase, meta->has("city") == false);

        CuAssertTrue(_testCase, meta->getMap().size() == 3);
    }
};

static void halMetaDataTest(CuTest *testCase) {
    MetaDataTest tester;
    tester.check(testCase);
}

static CuSuite *halMetaDataTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
// FIXME: commented this out till we implement this for mmap.
    if (false) {
        SUITE_ADD_TEST(suite, halMetaDataTest);
    }
    return suite;
}

int main(int argc, char *argv[]) {
    return runHalTestSuite(argc, argv, halMetaDataTestSuite());
}
