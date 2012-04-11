/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <string>
#include <iostream>
#include "halAlignmentTest.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignment.h"
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

TempCreateAlignment::TempCreateAlignment(AlignmentPtr alignment) :
  _alignment(alignment)
{
  _path = getTempFile();
  _alignment->createNew(_path);
}

TempCreateAlignment::~TempCreateAlignment()
{
  _alignment->close();
  removeTempFile(_path);
}

TempReadAlignment::TempReadAlignment(AlignmentPtr alignment, 
                                     char* path)
  : _path(path)
{
  alignment->open(_path, true);
  _alignment = alignment;
}

TempReadAlignment::~TempReadAlignment()
{
  _alignment->close();
}

void AlignmentTest::check(CuTest* testCase)
{
  _testCase = testCase;
  vector<AlignmentPtr> createInstances = getTestAlignmentInstances();
  vector<AlignmentPtr> readInstances = getTestAlignmentInstances();

  for (size_t i = 0; i < createInstances.size(); ++i)
  {
    TempCreateAlignment creater(createInstances[i]);
    createCallBack(creater._alignment);
    creater._alignment->close();
    
    TempReadAlignment checker(readInstances[i], creater._path);
    checkCallBack(checker._alignment);
  }
}

void AlignmentTestTrees::createCallBack(hal::AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  AlignmentTest::assertTrue(alignmentSize == 0);
}

void AlignmentTestTrees::checkCallBack(hal::AlignmentConstPtr alignment)
{

}

void halAlignmentTestTrees(CuTest *testCase)
{
  AlignmentTestTrees tester;
  tester.check(testCase);
}


CuSuite* halAlignmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halAlignmentTestTrees);
  return suite;
}

