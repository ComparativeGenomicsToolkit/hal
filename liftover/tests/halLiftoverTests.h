/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALLIFTOVERTESTS_H
#define _HALLIFTOVERTESTS_H

#include "hal.h"
#include "halApiTestSupport.h"

extern "C" {
#include "CuTest.h"
}

using namespace hal;

struct BedLiftoverTest : public AlignmentTest {
    private:
    void liftAndCheck(const Alignment *alignment,
                      const Genome *srcGenome,
                      const Genome *tgtGenome,
                      const std::string& inBed,
                      const std::string& expectBed,
                      bool outPSL = false, bool outPSLWithName = false);
    public:
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void testOneBranchLifts(const Alignment *alignment);
    void testMultiBranchLifts(const Alignment *alignment);
};

struct WiggleLiftoverTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void testOneBranchLifts(const Alignment *alignment);
    void testMultiBranchLifts(const Alignment *alignment);
};

CuSuite *halLiftoverTestSuite();

#endif
// Local Variables:
// mode: c++
// End:
