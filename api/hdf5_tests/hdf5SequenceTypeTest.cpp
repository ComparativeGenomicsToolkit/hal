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
#include <cstdlib>
#include <cassert>
#include <H5Cpp.h>
#include "allTests.h"
#include "hdf5ExternalArray.h"
#include "hdf5Sequence.h"
#include "hdf5BottomSegment.h"
#include "hdf5Test.h"
extern "C" {
#include "commonC.h"
}

using namespace H5;
using namespace hal;
using namespace std;

static const hsize_t chunkSizes[] = {0, 4, 8, 16, 32, 128, 512, 1000};
static const hsize_t numSizes = 8;
static const hsize_t maxNameLength[] = {1,10,15,123};
static const hsize_t numLengths = 4;

static void teardown()
{
  hdf5TestTeardown();
}

static void setup()
{
  hdf5TestSetup();
}

static string genName(hsize_t i, hsize_t maxLength)
{
  if (maxLength > 2 && i % 2 == 0)
  {
    maxLength /= 2;
  }
  string name;
  for (hsize_t i = 0; i < maxLength; ++i)
  {
    name.push_back((char)(48 + i % 50));
  }
  return name;
}

void hdf5SequenceTypeTest(CuTest *testCase)
{
  for (hsize_t lengthIdx = 0; lengthIdx < numLengths; ++lengthIdx)
  {
    hsize_t length = maxNameLength[lengthIdx];
    for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
    {
      hsize_t chunkSize = chunkSizes[chunkIdx];
      setup();
      try 
      {
        CompType datatype = HDF5Sequence::dataType(length);
        H5File file(H5std_string(fileName), H5F_ACC_TRUNC);

        HDF5ExternalArray myArray;
        DSetCreatPropList cparms;
        if (chunkSize > 0)
        {
          cparms.setChunk(1, &chunkSize);
        }
        myArray.create(&file, datasetName, datatype, N, cparms);
        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5Sequence sequence(NULL, &myArray, i);
          Sequence::Info seqInfo(genName(i, length), i * 2, i * 3, i * 4);
          sequence.set(i, seqInfo);
        }
        myArray.write();
        file.flush(H5F_SCOPE_LOCAL);
        file.close();

        H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
        HDF5ExternalArray readArray;
        readArray.load(&rfile, datasetName);

        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5Sequence sequence(NULL, &readArray, i);
          CuAssertTrue(testCase,
                       sequence.getName() == genName(i, length));
          CuAssertTrue(testCase, 
                       sequence.getStartPosition() == i);
          CuAssertTrue(testCase,
                       sequence.getSequenceLength() == i * 2);
          CuAssertTrue(testCase,
                       sequence.getNumTopSegments() == i * 3);
          CuAssertTrue(testCase,
                       sequence.getNumBottomSegments() == i * 4);
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

CuSuite* hdf5SequenceTypeTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, hdf5SequenceTypeTest);
  return suite;
}
