/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <string>
#include <iostream>
#include <cstdlib>
#include "halAlignmentTest.h"
#include "halAlignmentInstanceTest.h"
#include "halAlignment.h"
#include "halGenome.h"
extern "C" {
#include "commonC.h"
}

using namespace std;
using namespace hal;

string AlignmentTest::randomString(hal_size_t length)
{
  string s;
  s.resize(length);
  for (hal_size_t i = 0; i < length; ++i)
  {
    int r = rand() % 10;
    char c;
    switch (r) 
    {
    case 0 : c = 'a'; break;
    case 1 : c = 'c'; break;
    case 2 : c = 'g'; break;
    case 3 : c = 't'; break;
    case 4 : c = 'A'; break;
    case 5 : c = 'C'; break;
    case 6 : c = 'G'; break;
    case 7 : c = 'T'; break;
    case 8 : c = 'N'; break;
    case 9 : c = 'n'; break;
    default: c = '?'; break;
    }
    s[i] = c;
  }
  return s;
}

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
  alignment->close();
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
    _createPath = creater._path;
    createCallBack(creater._alignment);
    creater._alignment->close();
    
    TempReadAlignment checker(readInstances[i], creater._path);
    _checkPath = checker._path;
    checkCallBack(checker._alignment);
  }
}

void AlignmentTestTrees::createCallBack(hal::AlignmentPtr alignment)
{
  hal_size_t alignmentSize = alignment->getNumGenomes();
  CuAssertTrue(_testCase, alignmentSize == 0);
  
  alignment->addRootGenome("Root", 0);
  alignment->addLeafGenome("Leaf", "Root", 10);
  alignment->addRootGenome("NewRoot", 15);
  alignment->addLeafGenome("Leaf1", "Root", 4.1);
  alignment->addLeafGenome("Leaf2", "Root", 5.1);
  alignment->addLeafGenome("Leaf3", "Root", 6.1);
  alignment->addLeafGenome("Leaf4", "Root", 7.1);
  alignment->updateBranchLength("Root", "Leaf1", 3.0);
  alignment->updateBranchLength("Root", "Leaf2", 6.1);
  alignment->updateBranchLength("Root", "Leaf2", 5.1);
}

void AlignmentTestTrees::checkCallBack(hal::AlignmentConstPtr alignment)
{
  CuAssertTrue(_testCase, alignment->getRootName() == "NewRoot");
  CuAssertTrue(_testCase, alignment->getNewickTree() == 
        "((Leaf:10,Leaf1:3,Leaf2:5.1,Leaf3:6.1,Leaf4:7.1)Root:15)NewRoot;");
  CuAssertTrue(_testCase, alignment->getBranchLength("Root", "Leaf") == 10.0);
  CuAssertTrue(_testCase, alignment->getBranchLength("Root", "Leaf1") == 3.0);
  vector<string> children = alignment->getChildNames("Root");
  CuAssertTrue(_testCase, children.size() == 5);

  vector<string> leaves = alignment->getLeafNamesBelow("Leaf");
  CuAssertTrue(_testCase, leaves.size() == 0);
  leaves = alignment->getLeafNamesBelow("NewRoot");
  CuAssertTrue(_testCase, leaves.size() == 5);
  for (size_t i = 0; i < leaves.size(); ++i)
  {
    CuAssertTrue(_testCase, leaves[i][0] == 'L');
  }
  leaves = alignment->getLeafNamesBelow("Root");
  CuAssertTrue(_testCase, leaves.size() == 5);
  for (size_t i = 0; i < leaves.size(); ++i)
  {
    CuAssertTrue(_testCase, leaves[i][0] == 'L');
  }
}

void AlignmentTestEmpty::createCallBack(hal::AlignmentPtr alignment)
{

}

void AlignmentTestEmpty::checkCallBack(hal::AlignmentConstPtr alignment)
{
  CuAssertTrue(_testCase, alignment->getNumGenomes() == 0);
}

void AlignmentTestBadPathError::createCallBack(hal::AlignmentPtr alignment)
{
  
}

void AlignmentTestBadPathError::checkCallBack(hal::AlignmentConstPtr alignment)
{
#ifndef ENABLE_UDC
  // note: udc just exits here so test wont work
  string badPath(_checkPath);
  badPath += "blarg";
  try {
    ((Alignment*)alignment.get())->open(badPath, true);
  }
  catch(...)
  {
    return;
  }
  CuAssertTrue(_testCase, false);
#endif
}

void halAlignmentTestTrees(CuTest *testCase)
{
  AlignmentTestTrees tester;
  tester.check(testCase);
}

void halAlignmentTestEmpty(CuTest *testCase)
{
  AlignmentTestEmpty tester;
  tester.check(testCase);
}

void halAlignmentTestBadPathError(CuTest *testCase)
{
  AlignmentTestBadPathError tester;
  tester.check(testCase);
}

CuSuite* halAlignmentTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, halAlignmentTestTrees);
  SUITE_ADD_TEST(suite, halAlignmentTestEmpty);
  SUITE_ADD_TEST(suite, halAlignmentTestBadPathError);
  return suite;
}

