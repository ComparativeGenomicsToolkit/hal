/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "allTests.h"
#include "halAlignmentInstance.h"
#include "halAlignmentTest.h"
#include <cstdio>
#include <iostream>

int halRunAllTests(void) {
    CuString *output = CuStringNew();
    CuSuite *suite = CuSuiteNew();
    CuSuiteAddSuite(suite, halAlignmentTestSuite());
    CuSuiteAddSuite(suite, halMetaDataTestSuite());
    CuSuiteAddSuite(suite, halGenomeTestSuite());
    CuSuiteAddSuite(suite, halTopSegmentTestSuite());
    CuSuiteAddSuite(suite, halBottomSegmentTestSuite());
    CuSuiteAddSuite(suite, halSequenceTestSuite());
    CuSuiteAddSuite(suite, halColumnIteratorTestSuite());
    CuSuiteAddSuite(suite, halGappedSegmentIteratorTestSuite());
    CuSuiteAddSuite(suite, halRearrangementTestSuite());
    CuSuiteAddSuite(suite, halMappedSegmentTestSuite());
    CuSuiteAddSuite(suite, halValidateTestSuite());
    CuSuiteRun(suite);
    CuSuiteSummary(suite, output);
    CuSuiteDetails(suite, output);
    printf("%s\n", output->buffer);
    return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
    if (argc > 2) {
        std::cerr << "wrong # args: " << argv[0] << " [storageDriver]" << std::endl;
        ::exit(1);
    } else if (argc == 2) {
        storageDriverToTest = argv[1];
        if (not((storageDriverToTest == hal::STORAGE_FORMAT_HDF5) or (storageDriverToTest == hal::STORAGE_FORMAT_MMAP))) {
            std::cerr << "Invalid storage driver '" << storageDriverToTest << "', expected on of: " << hal::STORAGE_FORMAT_HDF5
                      << " or " << hal::STORAGE_FORMAT_MMAP << std::endl;
            ::exit(1);
        }
    }
    return halRunAllTests();
}
