/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCHAINGETBLOCKSTEST_H
#define _HALCHAINGETBLOCKSTEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "hal.h"

typedef std::set<hal::MappedSegmentConstPtr,
                 hal::MappedSegment::LessSource> MSRefSet;

struct ChainGetBlocksSimpleTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   virtual void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksInversionTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksOffsetTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksInversionOffsetTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksOffsetQRefTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksInversionOffsetQRefTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksInversionOffsetQSisTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ChainGetBlocksSimpleLiftoverTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(hal::AlignmentConstPtr alignment);
};



#endif
