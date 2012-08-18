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

void addIdenticalParentChild(hal::AlignmentPtr alignment,
                             size_t numSequences,
                             size_t numSegmentsPerSequence,
                             size_t segmentLength);

void makeInsertion(hal::BottomSegmentIteratorPtr bi);

void makeInsGap(hal::TopSegmentIteratorPtr ti);

void makeDelGap(hal::BottomSegmentIteratorPtr bi);

void makeInversion(hal::TopSegmentIteratorPtr ti, hal_size_t len);

struct RearrangementInsertionTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct RearrangementSimpleInversionTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct RearrangementGappedInversionTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};


#endif
