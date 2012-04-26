/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALVALIDATETEST_H
#define _HALVALIDATETEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

struct ValidateSmallTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ValidateMediumTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ValidateLargeTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

#endif
