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

using namespace hal;

struct GenomeMetaTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct GenomeCreateTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct GenomeUpdateTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
};

struct GenomeStringTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   std::string _string;
};

struct GenomeCopyTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   std::string _path;
   AlignmentPtr _secondAlignment;
};

struct GenomeCopySegmentsWhenSequencesOutOfOrderTest : public AlignmentTest
{
   void createCallBack(Alignment* alignment);
   void checkCallBack(const Alignment* alignment);
   std::string _path;
   AlignmentPtr _secondAlignment;
};

#endif
// Local Variables:
// mode: c++
// End:
