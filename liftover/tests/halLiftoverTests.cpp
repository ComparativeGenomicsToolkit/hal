/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cstdio>
#include "halLiftoverTests.h"



CuSuite* halLiftoverTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
//  SUITE_ADD_TEST(suite, halTopSegmentSimpleIteratorTest);
  return suite;
}

int halLiftoverRunAllTests(void) {
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();
  CuSuiteAddSuite(suite, halLiftoverTestSuite());
  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("%s\n", output->buffer);
  return suite->failCount > 0;
}

int main(int argc, char *argv[]) {
   
  return halLiftoverRunAllTests();
}
