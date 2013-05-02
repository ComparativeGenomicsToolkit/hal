/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAPPEDSEGMENTTEST_H
#define _HALMAPPEDSEGMENTTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

struct MappedSegmentMapUpTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   virtual void checkCallBack(hal::AlignmentConstPtr alignment);
   void testTopSegment(hal::AlignmentConstPtr alignment,
                       hal::TopSegmentIteratorConstPtr top);
};

struct MappedSegmentMapDownTest : public MappedSegmentMapUpTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void testBottomSegment(hal::AlignmentConstPtr alignment,
                          hal::BottomSegmentIteratorConstPtr bottom,
                          hal_size_t childIndex);

};

struct MappedSegmentMapAcrossTest : public MappedSegmentMapUpTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void testTopSegment(hal::AlignmentConstPtr alignment,
                       hal::TopSegmentIteratorConstPtr top);
};

struct MappedSegmentMapDupeTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};


#endif
