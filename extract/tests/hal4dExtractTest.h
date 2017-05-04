#ifndef _HAL4DEXTRACTTEST_H
#define _HAL4DEXTRACTTEST_H
#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

struct Bed4dExtractTest : public AlignmentTest
{
  void createCallBack(hal::AlignmentPtr alignment);
  void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ConservedBed4dExtractTest : public AlignmentTest
{
  void createCallBack(hal::AlignmentPtr alignment);
  void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct Block4dExtractTest : public AlignmentTest
{
  void createCallBack(hal::AlignmentPtr alignment);
  void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct ConservedBlock4dExtractTest : public AlignmentTest
{
  void createCallBack(hal::AlignmentPtr alignment);
  void checkCallBack(hal::AlignmentConstPtr alignment);
};

struct CDS4dExtractTest : public AlignmentTest
{
  void createCallBack(hal::AlignmentPtr alignment);
  void checkCallBack(hal::AlignmentConstPtr alignment);
};

#endif
