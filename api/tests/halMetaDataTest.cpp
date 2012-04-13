/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <string>
#include <iostream>
#include "halMetaDataTest.h"
#include "halAlignmentTest.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignment.h"
#include "halGenome.h"
#include "halMetaData.h"
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

void MetaDataTest::createCallBack(hal::AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  MetaData* meta = alignment->getMetaData();
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

void MetaDataTest::checkCallBack(hal::AlignmentConstPtr alignment)
{
  const MetaData* meta = alignment->getMetaData();
  
  CuAssertTrue(_testCase, meta->get("colour") == "black");
  CuAssertTrue(_testCase, meta->get("number") == "1");
  CuAssertTrue(_testCase, meta->get("animal") == "cat");
  
  CuAssertTrue(_testCase, meta->has("colour") == true);
  CuAssertTrue(_testCase, meta->has("city") == false);
  
  CuAssertTrue(_testCase, meta->getMap().size() == 3);
}

void halMetaDataTest(CuTest *testCase)
{
  MetaDataTest tester;
  tester.check(testCase);
}

CuSuite* halMetaDataTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halMetaDataTest);
  return suite;
}

