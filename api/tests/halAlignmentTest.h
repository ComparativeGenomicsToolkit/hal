/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALALIGNMENTTEST_H
#define _HALALIGNMENTTEST_H

#include <vector>
#include "halAlignment.h"
#include "allTests.h"

struct TempCreateAlignment 
{
   TempCreateAlignment(hal::AlignmentPtr alignment);
   char* _path;
   hal::AlignmentPtr _alignment;
   ~TempCreateAlignment();
};

struct TempReadAlignment 
{
   TempReadAlignment(hal::AlignmentPtr alignment,
                     char* path);
   char* _path;
   hal::AlignmentPtr _alignment;
   ~TempReadAlignment();
};

struct AlignmentTest
{
   AlignmentTest(){}
   virtual ~AlignmentTest(){}
   void check(CuTest *testCase);   
   virtual void createCallBack(hal::AlignmentPtr alignment) {}
   virtual void checkCallBack(hal::AlignmentConstPtr alignment) {}
   CuTest* _testCase;
   void assertTrue(bool b) {
     CuAssertTrue(_testCase, b);
   }
};

struct AlignmentTestTrees : public AlignmentTest
{
   AlignmentTestTrees() : AlignmentTest(){}
   ~AlignmentTestTrees(){}
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};
#endif
