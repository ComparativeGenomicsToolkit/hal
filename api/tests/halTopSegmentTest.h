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
   hal_size_t _length;
   hal_index_t _startPosition;
   hal_index_t _nextParalogyIndex;
   hal_index_t _parentIndex;
   bool _parentReversed;
   hal_index_t _arrayIndex;
   hal_index_t _bottomParseIndex;
   void setRandom();
   void set(hal_index_t startPosition,
            hal_size_t length,
            hal_index_t parentIndex = hal::NULL_INDEX,
            bool parentReversed = false,
            hal_index_t bottomParseIndex = hal::NULL_INDEX,
            hal_index_t nextParalogyIndex = hal::NULL_INDEX);
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

struct TopSegmentIteratorParseTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct TopSegmentIteratorToSiteTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   void checkGenome(const hal::Genome* genome);
};

struct TopSegmentIteratorReverseTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct TopSegmentIsGapTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

#endif
