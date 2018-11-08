#ifndef _HAL4DEXTRACTTEST_H
#define _HAL4DEXTRACTTEST_H
#include "halAlignmentTest.h"

extern "C" {
#include "CuTest.h"
}

using namespace hal;

struct Bed4dExtractTest : public AlignmentTest
{
  void createCallBack(Alignment* alignment);
  void checkCallBack(const Alignment* alignment);
};

struct ConservedBed4dExtractTest : public AlignmentTest
{
  void createCallBack(Alignment* alignment);
  void checkCallBack(const Alignment* alignment);
};

struct Block4dExtractTest : public AlignmentTest
{
  void createCallBack(Alignment* alignment);
  void checkCallBack(const Alignment* alignment);
};

struct ConservedBlock4dExtractTest : public AlignmentTest
{
  void createCallBack(Alignment* alignment);
  void checkCallBack(const Alignment* alignment);
};

struct CDS4dExtractTest : public AlignmentTest
{
  void createCallBack(Alignment* alignment);
  void checkCallBack(const Alignment* alignment);
};

#endif
