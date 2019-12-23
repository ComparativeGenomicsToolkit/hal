/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFTESTS_H
#define _HALMAFTESTS_H

#include "halApiTestSupport.h"

extern "C" {
#include "CuTest.h"
}

CuSuite *halMafExportTestSuite();
CuSuite *halMafBlockTestSuite();

#endif
// Local Variables:
// mode: c++
// End:
