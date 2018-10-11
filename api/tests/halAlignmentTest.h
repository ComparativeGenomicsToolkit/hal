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

class AlignmentTest
{
public:
   AlignmentTest() {
   }
   virtual ~AlignmentTest() {
   }
   void check(CuTest *testCase);   
   virtual void createCallBack(hal::AlignmentPtr alignment) {}
   virtual void checkCallBack(hal::AlignmentConstPtr alignment) {}
   CuTest* _testCase;
   std::string _createPath;
   std::string _checkPath;
   static std::string randomString(hal_size_t length);
   void checkOne(CuTest* testCase,
                 const std::string& storageFormat);
};

#endif
