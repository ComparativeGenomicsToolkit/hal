/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALREARRANGEMENTTEST_H
#define _HALREARRANGEMENTTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"
#include "allTests.h"

using namespace hal;

void addIdenticalParentChild(Alignment* alignment,
                             size_t numSequences,
                             size_t numSegmentsPerSequence,
                             size_t segmentLength);

void makeInsertion(BottomSegmentIteratorPtr bi);

void makeInsGap(TopSegmentIteratorPtr ti);

void makeDelGap(BottomSegmentIteratorPtr bi);

void makeInversion(TopSegmentIteratorPtr ti, hal_size_t len);

struct RearrangementInsertionTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct RearrangementSimpleInversionTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct RearrangementGappedInversionTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};


#endif
