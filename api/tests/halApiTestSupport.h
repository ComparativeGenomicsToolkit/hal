/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALAPITESTSUPPORT_H
#define _HALAPITESTSUPPORT_H

#include "halAlignmentInstance.h"
#include <vector>
extern "C" {
#include "CuTest.h"
#include "commonC.h"
}

using namespace hal;
using namespace std;

Alignment *getTestAlignmentInstances(const string &storageFormat, const string &alignmentPath, unsigned mode);

/** parse command line and run a test suite for the given storage driver,
 * return exit code  */
int runHalTestSuite(int argc, char *argv[], CuSuite *suite);

/* Base class for alignment tests.  Handles setup of test HAL and has required
 * methods. */
class AlignmentTest {
  public:
    AlignmentTest() {
    }
    virtual ~AlignmentTest() {
    }
    void check(CuTest *testCase);
    virtual void createCallBack(Alignment *alignment) {
    }
    virtual void checkCallBack(const Alignment *alignment) {
    }
    CuTest *_testCase;
    string _createPath;
    string _checkPath;
    static string randomString(hal_size_t length);
    void checkOne(CuTest *testCase, const string &storageFormat);
};



#endif
// Local Variables:
// mode: c++
// End:
