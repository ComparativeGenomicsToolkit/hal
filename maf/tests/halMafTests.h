/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFTESTS_H
#define _HALMAFTESTS_H

#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

CuSuite *halMafExportTestSuite();
CuSuite *halMafBlockTestSuite();

#endif
