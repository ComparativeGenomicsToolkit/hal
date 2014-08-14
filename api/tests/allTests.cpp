/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include "allTests.h"

int halRunAllTests(void) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew(); 
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
   
  return halRunAllTests();
}
