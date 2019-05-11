/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOLUMNITERATORTEST_H
#define _HALCOLUMNITERATORTEST_H

#include "allTests.h"
#include "hal.h"
#include "halAlignmentTest.h"
#include <vector>

using namespace hal;

struct ColumnIteratorBaseTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct ColumnIteratorDepthTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void checkGenome(const Genome *genome);
};

struct ColumnIteratorDupTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void checkGenome(const Genome *genome);
};

struct ColumnIteratorInvTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void checkGenome(const Genome *genome);
};

struct ColumnIteratorGapTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct ColumnIteratorMultiGapTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct ColumnIteratorMultiGapInvTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct ColumnIteratorPositionCacheTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

#endif
// Local Variables:
// mode: c++
// End:
