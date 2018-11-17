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

using namespace hal;

struct MappedSegmentMapUpTest : virtual public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   virtual void checkCallBack(const Alignment* alignment);
   void testTopSegment(const Alignment* alignment,
                       TopSegmentIteratorPtr top,
                       const std::string& ancName);
};

struct MappedSegmentParseTest : virtual public MappedSegmentMapUpTest
{
   virtual void checkCallBack(const Alignment* alignment);
   void testTopSegment(const Alignment* alignment, 
                       TopSegmentIteratorPtr top);                   
};

struct MappedSegmentMapDownTest : public MappedSegmentMapUpTest
{
   void checkCallBack(const Alignment* alignment);
   void testBottomSegment(const Alignment* alignment,
                          BottomSegmentIteratorPtr bottom,
                          hal_size_t childIndex);

};

struct MappedSegmentMapAcrossTest : public MappedSegmentMapUpTest
{
   void checkCallBack(const Alignment* alignment);
   void testTopSegment(const Alignment* alignment,
                       TopSegmentIteratorPtr top);
};

struct MappedSegmentMapDupeTest : virtual public AlignmentTest
{
   virtual void createCallBack(Alignment* alignment);
   virtual void checkCallBack(const Alignment* alignment);
};

struct MappedSegmentMapExtraParalogsTest : virtual public AlignmentTest
{
  virtual void createCallBack(Alignment* alignment);
  virtual void checkCallBack(const Alignment* alignment);
};

struct MappedSegmentColCompareTest : virtual public AlignmentTest
{
   virtual void createCallBack(Alignment* alignment) = 0;
   void checkCallBack(const Alignment* alignment);
   void createColArray();
   void createBlockArray();
   void compareArrays();
   std::vector<std::set<std::pair<hal_index_t, bool> > >_colArray;
   std::vector<std::set<std::pair<hal_index_t, bool> > >_blockArray;
   const Genome* _ref;
   const Genome* _tgt;
};

struct 
MappedSegmentColCompareTestCheck1 : virtual public MappedSegmentMapUpTest,
                                    virtual public MappedSegmentColCompareTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct 
MappedSegmentColCompareTestCheck2 : virtual public MappedSegmentMapDupeTest,
                                    virtual public MappedSegmentColCompareTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct MappedSegmentColCompareTest1 : public MappedSegmentColCompareTest
{
   void createCallBack(Alignment* alignment);
};

struct MappedSegmentColCompareTest2 : public MappedSegmentColCompareTest
{
   void createCallBack(Alignment* alignment);
};

struct MappedSegmentColCompareTest3 : public MappedSegmentColCompareTest
{
   void createCallBack(Alignment* alignment);
};


#endif
// Local Variables:
// mode: c++
// End:
