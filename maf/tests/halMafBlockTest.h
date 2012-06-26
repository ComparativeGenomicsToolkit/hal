
/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFBLOCKTEST_H
#define _HALMAFBLOCKTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "halMafTests.h"

struct MafBlockCreateTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};


#endif
