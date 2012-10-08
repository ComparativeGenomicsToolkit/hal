/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include <H5Cpp.h>
#include "allTests.h"

int hdf5RunAllTests(void) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();
  CuSuiteAddSuite(suite, hdf5TestSuite());
  CuSuiteAddSuite(suite, hdf5ExternalArrayTestSuite());
  CuSuiteAddSuite(suite, hdf5DNATypeTestSuite());
  CuSuiteAddSuite(suite, hdf5SegmentTypeTestSuite());
  CuSuiteAddSuite(suite, hdf5SequenceTypeTestSuite());
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
  return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
  H5::Exception::dontPrint();
  return hdf5RunAllTests();
}
