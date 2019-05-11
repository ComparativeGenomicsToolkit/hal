/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALBOTTOMSEGMENTTEST_H
#define _HALBOTTOMSEGMENTTEST_H

#include "allTests.h"
#include "hal.h"
#include "halAlignmentTest.h"
#include <vector>

using namespace hal;

struct BottomSegmentStruct {
    hal_size_t _length;
    hal_index_t _startPosition;
    std::vector<std::pair<hal_index_t, bool>> _children;
    hal_index_t _arrayIndex;
    hal_index_t _topParseIndex;
    void setRandom(hal_size_t numChildren);
    void applyTo(BottomSegmentIteratorPtr it) const;
    void compareTo(BottomSegmentIteratorPtr it, CuTest *testCase) const;
    void set(hal_index_t startPosition, hal_size_t length, hal_index_t topParseIndex = NULL_INDEX);
};

struct BottomSegmentSimpleIteratorTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    std::vector<BottomSegmentStruct> _bottomSegments;
};

struct BottomSegmentSequenceTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    std::vector<BottomSegmentStruct> _bottomSegments;
    std::vector<std::string> _sequences;
};

struct BottomSegmentIteratorParseTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct BottomSegmentIteratorToSiteTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
    void checkGenome(const Genome *genome);
};

struct BottomSegmentIteratorReverseTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

struct BottomSegmentIsGapTest : public AlignmentTest {
    void createCallBack(Alignment *alignment);
    void checkCallBack(const Alignment *alignment);
};

#endif
// Local Variables:
// mode: c++
// End:
