
/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALSEQUENCETEST_H
#define _HALSEQUENCETEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

using namespace hal;

struct SequenceCreateTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct SequenceIteratorTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct SequenceUpdateTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

#endif
// Local Variables:
// mode: c++
// End:
