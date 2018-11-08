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

using namespace hal;

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
            hal_index_t parentIndex = NULL_INDEX,
            bool parentReversed = false,
            hal_index_t bottomParseIndex = NULL_INDEX,
            hal_index_t nextParalogyIndex = NULL_INDEX);
   void applyTo(TopSegmentIteratorPtr it) const;
   void compareTo(TopSegmentIteratorPtr it,
     CuTest* testCase) const;
};

struct TopSegmentSimpleIteratorTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   std::vector<TopSegmentStruct> _topSegments;
};

struct TopSegmentSequenceTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   std::vector<TopSegmentStruct> _topSegments;
   std::vector<std::string> _sequences;
};

struct TopSegmentIteratorParseTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct TopSegmentIteratorToSiteTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   void checkGenome(const Genome* genome);
};

struct TopSegmentIteratorReverseTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct TopSegmentIsGapTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

#endif
