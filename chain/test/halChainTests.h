/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCHAINTESTS_H
#define _HALCHAINTESTS_H

#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

CuSuite *halChainGetBlocksTestSuite();

#endif
