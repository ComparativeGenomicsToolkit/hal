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

using namespace hal;

typedef std::set<MappedSegmentPtr,
                 MappedSegment::LessSource> MSRefSet;

struct ChainGetBlocksSimpleTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   virtual void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksInversionTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksOffsetTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksInversionOffsetTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksOffsetQRefTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksInversionOffsetQRefTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksInversionOffsetQSisTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};

struct ChainGetBlocksSimpleLiftoverTest : public ChainGetBlocksSimpleTest
{
   void checkCallBack(const Alignment* alignment);
};



#endif
