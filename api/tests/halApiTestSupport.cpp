/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halApiTestSupport.h"
#include "halAlignmentInstance.h"
#include <H5Cpp.h>
#include <iostream>
#include <cstdlib>
#include <unistd.h>


using namespace std;
using namespace hal;

/* Global set from command line containing the storage drives to use.  * This
 * allows parallelizing tests.
 */
static string storageDriverToTest;

Alignment *getTestAlignmentInstances(const std::string &storageFormat, const std::string &alignmentPath, unsigned mode) {
    if (storageFormat == STORAGE_FORMAT_HDF5) {
        return hdf5AlignmentInstance(alignmentPath, mode, hdf5DefaultFileCreatPropList(), hdf5DefaultFileAccPropList(),
                                     hdf5DefaultDSetCreatPropList());

    } else if (storageFormat == hal::STORAGE_FORMAT_MMAP) {
        // We use a default init size of only 1GiB here, because the test
        // alignments we create are relatively small.
        return mmapAlignmentInstance(alignmentPath, mode, 1024 * 1024 * 1024);
    } else {
        throw hal_exception("invalid storage format: " + storageFormat);
    }
}

/** parse command line and run a test suite for the given storage driver,
 * return exit code  */
int runHalTestSuite(int argc, char *argv[], CuSuite *suite) {
    if (argc > 2) {
        cerr << "wrong # args: " << argv[0] << " [storageDriver]" << endl;
        return 1;
    } else if (argc == 2) {
        storageDriverToTest = argv[1];
        if (not((storageDriverToTest == hal::STORAGE_FORMAT_HDF5) or (storageDriverToTest == hal::STORAGE_FORMAT_MMAP))) {
            cerr << "Invalid storage driver '" << storageDriverToTest << "', expected on of: " << hal::STORAGE_FORMAT_HDF5
                      << " or " << hal::STORAGE_FORMAT_MMAP << endl;
            return 1;
        }
    } else {
        storageDriverToTest = hal::STORAGE_FORMAT_HDF5;
    }

    CuString *output = CuStringNew();
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    cerr << argv[0] << output->buffer << endl;
    return (suite->failCount > 0) ? 1 : 0;

}

string AlignmentTest::randomString(hal_size_t length) {
    string s;
    s.resize(length);
    for (hal_size_t i = 0; i < length; ++i) {
        int r = rand() % 10;
        char c;
        switch (r) {
        case 0:
            c = 'a';
            break;
        case 1:
            c = 'c';
            break;
        case 2:
            c = 'g';
            break;
        case 3:
            c = 't';
            break;
        case 4:
            c = 'A';
            break;
        case 5:
            c = 'C';
            break;
        case 6:
            c = 'G';
            break;
        case 7:
            c = 'T';
            break;
        case 8:
            c = 'N';
            break;
        case 9:
            c = 'n';
            break;
        default:
            c = '?';
            break;
        }
        s[i] = c;
    }
    return s;
}

void AlignmentTest::check(CuTest *testCase) {
    _testCase = testCase;
    try {
        if (storageDriverToTest.empty() or (storageDriverToTest == STORAGE_FORMAT_HDF5)) {
            checkOne(testCase, STORAGE_FORMAT_HDF5);
        }
        if (storageDriverToTest.empty() or (storageDriverToTest == STORAGE_FORMAT_MMAP)) {
            checkOne(testCase, STORAGE_FORMAT_MMAP);
        }
    } catch (const exception &e) {
        CuFail(testCase, stString_print("Caught exception while testing: %s", e.what()));
    }
}

void AlignmentTest::checkOne(CuTest *testCase, const string &storageFormat) {
    string alignmentPath = getTempFile();

    // test with created
    AlignmentPtr calignment(getTestAlignmentInstances(storageFormat, alignmentPath, CREATE_ACCESS));
    _createPath = alignmentPath;
    createCallBack(calignment.get());
    calignment->close();

    // test with existing alignment
    AlignmentPtr ralignment(getTestAlignmentInstances(storageFormat, alignmentPath, READ_ACCESS));
    _checkPath = alignmentPath;
    checkCallBack(ralignment.get());
    ralignment->close();

    ::unlink(alignmentPath.c_str());
}

