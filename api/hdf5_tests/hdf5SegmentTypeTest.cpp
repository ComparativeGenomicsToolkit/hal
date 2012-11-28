/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * Test the Segment offset nonsense for HDF5 arrays
 */

#include <iostream>
#include <string>
#include <cassert>
#include <H5Cpp.h>
#include "allTests.h"
#include "hdf5ExternalArray.h"
#include "hdf5TopSegment.h"
#include "hdf5BottomSegment.h"
#include "hdf5Test.h"
extern "C" {
#include "commonC.h"
}

using namespace H5;
using namespace hal;
using namespace std;

static const hsize_t chunkSizes[] = {0, 4, N/10, N/2, N};
static const hsize_t numSizes = 5;
static const hsize_t numChildrens[] = {0,1,2,3,50};
static const hsize_t numNumChildrens = 5;

static void teardown()
{
  hdf5TestTeardown();
}

static void setup()
{
  hdf5TestSetup();
}

void hdf5SegmentTypeTestTop(CuTest *testCase)
{
  for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
  {
    hsize_t chunkSize = chunkSizes[chunkIdx];
    setup();
    try 
    {
      CompType datatype = HDF5TopSegment::dataType();
      H5File file(H5std_string(fileName), H5F_ACC_TRUNC);

      HDF5ExternalArray myArray;
      DSetCreatPropList cparms;
      if (chunkSize > 0)
      {
        cparms.setChunk(1, &chunkSize);
      }
      myArray.create(&file, datasetName, datatype, N + 1, &cparms);
      for (hsize_t i = 0; i < N; ++i)
      {
        HDF5TopSegment segment(NULL, &myArray, i);
        segment.setCoordinates(i, 1);
        segment.setNextParalogyIndex(i * 3 + 1);
        segment.setParentIndex(i * 4);
        segment.setParentReversed(i % 2 ? true : false);
        segment.setBottomParseIndex(i * 5);
      }
      myArray.write();
      file.flush(H5F_SCOPE_LOCAL);
      file.close();

      H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
      HDF5ExternalArray readArray;
      readArray.load(&rfile, datasetName);
      for (hsize_t i = 0; i < N; ++i)
      {
         HDF5TopSegment segment(NULL, &readArray, i);
         CuAssertTrue(testCase, 
                      segment.getStartPosition() == (hal_index_t)i);
         CuAssertTrue(testCase,
                      segment.getLength() == 1);
         CuAssertTrue(testCase,
                      segment.getNextParalogyIndex() == (hal_index_t)i * 3 + 1);
         CuAssertTrue(testCase,
                      segment.getParentIndex() == (hal_index_t)i * 4);
         CuAssertTrue(testCase,
                      segment.getParentReversed() == i % 2 ? true : false);
         CuAssertTrue(testCase,
                      segment.getBottomParseIndex() == (hal_index_t)i * 5);
      }
    }
    catch(Exception& exception)
    {
      cerr << exception.getCDetailMsg() << endl;
      CuAssertTrue(testCase, 0);
    }
    catch(...)
    {
      CuAssertTrue(testCase, 0);
    }
    teardown();
  }
}

void hdf5SegmentTypeTestBottom(CuTest *testCase)
{
  for (hsize_t numChildIdx = 0; numChildIdx < numNumChildrens; ++numChildIdx)
  {
    hsize_t numChildren = numChildrens[numChildIdx];
    for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
    {
      hsize_t chunkSize = chunkSizes[chunkIdx];
      setup();
      try 
      {
        CompType datatype = HDF5BottomSegment::dataType(numChildren);
        H5File file(H5std_string(fileName), H5F_ACC_TRUNC);

        HDF5ExternalArray myArray;
        DSetCreatPropList cparms;
        if (chunkSize > 0)
        {
          cparms.setChunk(1, &chunkSize);
        }
        myArray.create(&file, datasetName, datatype, N + 1, &cparms);
        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5BottomSegment segment(NULL, &myArray, i);
          segment.setCoordinates(i, 1);
          segment.setTopParseIndex(i * 5);
          for (hsize_t j = 0; j < numChildren; ++j)
          {
            segment.setChildIndex(j, i * 4 + j);
            segment.setChildReversed(j, (i+j) % 2 ? true : false);
          }
        }
        myArray.write();
        file.flush(H5F_SCOPE_LOCAL);
        file.close();

        H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
        HDF5ExternalArray readArray;
        readArray.load(&rfile, datasetName);

        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5BottomSegment segment(NULL, &readArray, i);
          CuAssertTrue(testCase, 
                       segment.getStartPosition() == (hal_index_t)i * 1);
          CuAssertTrue(testCase,
                       segment.getLength() == 1);
          CuAssertTrue(testCase,
                       segment.getTopParseIndex() == (hal_index_t)i * 5);
          for (hsize_t j = 0; j < numChildren; ++j)
          {
            CuAssertTrue(testCase,
                         segment.getChildIndex(j) == (hal_index_t)(i * 4 + j));
            CuAssertTrue(testCase,
                         segment.getChildReversed(j) == (i+j) % 2 ? true : false);
          }

        }
      }
      catch(Exception& exception)
      {
        cerr << exception.getCDetailMsg() << endl;
        CuAssertTrue(testCase, 0);
      }
      catch(...)
      {
        CuAssertTrue(testCase, 0);
      }
      teardown();
    }
  }
}

CuSuite* hdf5SegmentTypeTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, hdf5SegmentTypeTestTop);
  SUITE_ADD_TEST(suite, hdf5SegmentTypeTestBottom);
  return suite;
}
