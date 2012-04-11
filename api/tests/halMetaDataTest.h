/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMETADATATEST_H
#define _HALMETADATATEST_H

#include <vector>
#include "halAlignmentTest.h"
#include "halMetaData.h"
#include "allTests.h"

struct MetaDataTest : public AlignmentTest
{
   void createCallBack(hal::AlignmentPtr alignment);
   void checkCallBack(hal::AlignmentConstPtr alignment);
};

#endif
