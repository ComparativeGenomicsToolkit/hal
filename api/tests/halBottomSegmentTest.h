/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENTTEST_H
#define _HALBOTTOMSEGMENTTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

struct BottomSegmentStruct {
   hal_index_t _length;
   hal_index_t _startPosition;
   hal_index_t _nextParalogyIndex;
   std::vector<std::pair<hal_index_t, hal_bool_t> >_children;
   hal_index_t _arrayIndex;
   hal_index_t _topParseIndex;
   hal_index_t _topParseOffset;
   void setRandom(hal_size_t numChildren);
   void applyTo(hal::BottomSegmentIteratorPtr it) const;
   void compareTo(hal::BottomSegmentIteratorConstPtr it,
     CuTest* testCase) const;
};

struct BottomSegmentSimpleIteratorTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::vector<BottomSegmentStruct> _bottomSegments;
};

struct BottomSegmentSequenceTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::vector<BottomSegmentStruct> _bottomSegments;
   std::vector<std::string> _sequences;
};

#endif
