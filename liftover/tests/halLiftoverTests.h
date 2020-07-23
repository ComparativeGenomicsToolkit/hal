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
    void liftAndCheck(AlignmentConstPtr alignment,
                      const Genome *srcGenome,
                      const Genome *tgtGenome,
                      const std::string& inBed,
                      const std::string& expectBed,
                      bool outPSL = false, bool outPSLWithName = false);
    public:
    void createCallBack(AlignmentPtr alignment);
    void checkCallBack(AlignmentConstPtr alignment);
    void testOneBranchLifts(AlignmentConstPtr alignment);
    void testMultiBranchLifts(AlignmentConstPtr alignment);
};

struct WiggleLiftoverTest : public AlignmentTest {
    void createCallBack(AlignmentPtr alignment);
    void checkCallBack(AlignmentConstPtr alignment);
    void testOneBranchLifts(AlignmentConstPtr alignment);
    void testMultiBranchLifts(AlignmentConstPtr alignment);
};

CuSuite *halLiftoverTestSuite();

#endif
// Local Variables:
// mode: c++
// End:
