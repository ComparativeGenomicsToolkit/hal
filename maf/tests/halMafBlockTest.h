
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

using namespace hal;

struct MafBlockCreateTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};


#endif
