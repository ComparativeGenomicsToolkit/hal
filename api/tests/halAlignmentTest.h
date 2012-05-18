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
   const char* _createPath;
   const char* _checkPath;
   static std::string randomString(hal_size_t length);
};

struct AlignmentTestTrees : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct AlignmentTestEmpty : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct AlignmentTestBadPathError : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};


#endif
