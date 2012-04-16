/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENTTEST_H
#define _HALTOPSEGMENTTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

struct TopSegmentStruct {
   hal_index_t _length;
   hal_index_t _startPosition;
   hal_index_t _nextParalogyIndex;
   hal_index_t _parentIndex;
   hal_index_t _arrayIndex;
   hal_index_t _bottomParseIndex;
   hal_index_t _bottomParseOffset;
   void setRandom();
   void applyTo(hal::TopSegmentIteratorPtr it) const;
   void compareTo(hal::TopSegmentIteratorConstPtr it,
     CuTest* testCase) const;
};

struct TopSegmentSimpleIteratorTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::vector<TopSegmentStruct> _topSegments;
};

struct TopSegmentSequenceTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::vector<TopSegmentStruct> _topSegments;
   std::vector<std::string> _sequences;
};

#endif
