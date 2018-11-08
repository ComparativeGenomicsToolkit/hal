/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDSEGMENTTEST_H
#define _HALGAPPEDSEGMENTTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "halTopSegmentTest.h"
#include "hal.h"
#include "allTests.h"

using namespace hal;

struct GappedSegmentSimpleIteratorTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct GappedSegmentSimpleIteratorTest2 : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct GappedSegmentIteratorIndelTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};



#endif
