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

struct MappedSegmentMapUpTest : virtual public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   virtual void checkCallBack(hal::AlignmentConstPtr alignment);
   void testTopSegment(hal::AlignmentConstPtr alignment,
                       hal::TopSegmentIteratorConstPtr top,
                       const std::string& ancName);
};

struct MappedSegmentParseTest : virtual public MappedSegmentMapUpTest
{
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

struct MappedSegmentMapDupeTest : virtual public AlignmentTest
{
   virtual void createCallBack(hal::AlignmentPtr alignment);
   virtual void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct MappedSegmentMapExtraParalogsTest : virtual public AlignmentTest
{
  virtual void createCallBack(hal::AlignmentPtr alignment);
  virtual void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct MappedSegmentColCompareTest : virtual public AlignmentTest
{
   virtual void createCallBack(hal::AlignmentPtr alignment) = 0;
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void createColArray();
   void createBlockArray();
   void compareArrays();
   std::vector<std::set<std::pair<hal_index_t, bool> > >_colArray;
   std::vector<std::set<std::pair<hal_index_t, bool> > >_blockArray;
   const hal::Genome* _ref;
   const hal::Genome* _tgt;
};

struct 
MappedSegmentColCompareTestCheck1 : virtual public MappedSegmentMapUpTest,
                                    virtual public MappedSegmentColCompareTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct 
MappedSegmentColCompareTestCheck2 : virtual public MappedSegmentMapDupeTest,
                                    virtual public MappedSegmentColCompareTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
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
