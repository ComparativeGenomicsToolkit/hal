#ifndef _HALGLOBALCOLUMNITERATORTEST_H
#define _HALGLOBALCOLUMNITERATORTEST_H
#include "allTests.h"
#include "halAlignmentTest.h"
#include "halAlignment.h"

struct GlobalColumnIteratorBasicTest : public AlignmentTest
{
    void createCallBack(hal::AlignmentPtr alignment);
    void checkCallBack(hal::AlignmentConstPtr alignment);
};

#endif
