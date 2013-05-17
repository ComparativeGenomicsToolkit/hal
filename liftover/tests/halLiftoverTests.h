/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLIFTOVERTESTS_H
#define _HALLIFTOVERTESTS_H

#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

CuSuite *halLiftoverTestSuite();

#endif
