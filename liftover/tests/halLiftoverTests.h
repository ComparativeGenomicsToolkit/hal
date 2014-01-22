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

struct BedLiftoverTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void testOneBranchLifts(hal::AlignmentConstPtr alignment);
   void testMultiBranchLifts(hal::AlignmentConstPtr alignment);
};

struct WiggleLiftoverTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void testOneBranchLifts(hal::AlignmentConstPtr alignment);
   void testMultiBranchLifts(hal::AlignmentConstPtr alignment);
};

CuSuite *halLiftoverTestSuite();

#endif
