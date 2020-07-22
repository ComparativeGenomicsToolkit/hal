
/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAFBLOCKTEST_H
#define _HALMAFBLOCKTEST_H

#include "hal.h"
#include "halApiTestSupport.h"
#include "halMafTests.h"
#include <vector>

using namespace hal;

struct MafBlockCreateTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment);
    void checkCallBack(AlignmentConstPtr alignment);
};

#endif
// Local Variables:
// mode: c++
// End:
