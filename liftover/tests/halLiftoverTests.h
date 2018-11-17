/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLIFTOVERTESTS_H
#define _HALLIFTOVERTESTS_H

#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

using namespace hal;

struct BedLiftoverTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   void testOneBranchLifts(const Alignment* alignment);
   void testMultiBranchLifts(const Alignment* alignment);
};

struct WiggleLiftoverTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   void testOneBranchLifts(const Alignment* alignment);
   void testMultiBranchLifts(const Alignment* alignment);
};

CuSuite *halLiftoverTestSuite();

#endif
// Local Variables:
// mode: c++
// End:
