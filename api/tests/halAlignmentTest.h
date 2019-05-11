/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTTEST_H
#define _HALALIGNMENTTEST_H

#include "allTests.h"
#include "halAlignment.h"
#include <string>
#include <vector>

using namespace hal;

/* Global set from command line containing the storage drives to use.  * This
 * allows parallelizing tests.  If empty, all storage drives are tested.
 */
extern std::string storageDriverToTest;

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
    std::string _createPath;
    std::string _checkPath;
    static std::string randomString(hal_size_t length);
    void checkOne(CuTest *testCase, const std::string &storageFormat);
};

#endif
// Local Variables:
// mode: c++
// End:
