/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGENOMETEST_H
#define _HALGENOMETEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "halGenome.h"
#include "allTests.h"

struct GenomeMetaTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct GenomeCreateTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct GenomeUpdateTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct GenomeStringTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::string _string;
};

struct GenomeCopyTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::string _path;
   hal::AlignmentPtr _secondAlignment;
};

struct GenomeCopySegmentsWhenSequencesOutOfOrderTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
   std::string _path;
   hal::AlignmentPtr _secondAlignment;
};

#endif
