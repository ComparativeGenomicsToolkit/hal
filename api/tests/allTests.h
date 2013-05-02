/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _ALLTESTS_H
#define _ALLTESTS_H

extern "C" {
#include "CuTest.h"
}

CuSuite *halAlignmentTestSuite();
CuSuite* halMetaDataTestSuite();
CuSuite* halGenomeTestSuite();
CuSuite* halTopSegmentTestSuite();
CuSuite* halBottomSegmentTestSuite();
CuSuite* halSequenceTestSuite();
CuSuite* halValidateTestSuite();
CuSuite* halColumnIteratorTestSuite();
CuSuite* halRearrangementTestSuite();
CuSuite* halMappedSegmentTestSuite();
CuSuite* halGappedSegmentIteratorTestSuite();

#endif
