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

struct MappedSegmentColCompareTest : public AlignmentTest
{
   virtual void createCallBack(hal::AlignmentPtr alignment) = 0;
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void createColArray(const hal::Genome* ref, const hal::Genome* tgt);
   void createBlockArray(const hal::Genome* ref, const hal::Genome* tgt);
   void compareArrays();
   std::vector<std::map<hal_index_t, bool> >_colArray;
   std::vector<std::map<hal_index_t, bool> >_blockArray;
};

struct MappedSegmentColCompareTest1 : public MappedSegmentColCompareTest
{
   void createCallBack(hal::AlignmentPtr alignment);
};

struct MappedSegmentColCompareTest2 : public MappedSegmentColCompareTest
{
   void createCallBack(hal::AlignmentPtr alignment);
};

struct MappedSegmentColCompareTest3 : public MappedSegmentColCompareTest
{
   void createCallBack(hal::AlignmentPtr alignment);
};


#endif
