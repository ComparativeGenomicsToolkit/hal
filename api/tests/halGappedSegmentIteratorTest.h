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

struct GappedSegmentSimpleIteratorTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct GappedSegmentSimpleIteratorTest2 : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct GappedSegmentIteratorIndelTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};



#endif
