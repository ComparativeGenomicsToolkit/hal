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

static const hsize_t chunkSizes[] = {0, 4, N/10, N/2, N};
static const hsize_t numSizes = 5;
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
        CompType idxDatatype = HDF5Sequence::idxDataType();
        StrType nameDatatype = HDF5Sequence::nameDataType(length);
        H5File file(H5std_string(fileName), H5F_ACC_TRUNC);

        HDF5ExternalArray myIdxArray;
        HDF5ExternalArray myNameArray;
        DSetCreatPropList cparms;
        if (chunkSize > 0)
        {
          cparms.setChunk(1, &chunkSize);
        }
        myIdxArray.create(&file, datasetName, idxDatatype, N + 1, &cparms);
        myNameArray.create(&file, datasetName + "_n", nameDatatype, N, 
                           &cparms);
        hal_size_t totalTopSegments = 0;
        hal_size_t totalBottomSegments = 0;
        hal_index_t startPosition = 0;
        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5Sequence sequence(NULL, &myIdxArray, &myNameArray, i);
          Sequence::Info seqInfo(genName(i, length), i * 2, i * 3, i * 4);
          sequence.set(startPosition, seqInfo, totalTopSegments, 
                       totalBottomSegments);
          startPosition += i * 2;
          totalTopSegments += seqInfo._numTopSegments;
          totalBottomSegments += seqInfo._numBottomSegments;
        }
        myIdxArray.write();
        myNameArray.write();
        file.flush(H5F_SCOPE_LOCAL);
        file.close();

        H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
        HDF5ExternalArray readIdxArray;
        HDF5ExternalArray readNameArray;
        readIdxArray.load(&rfile, datasetName);
        readNameArray.load(&rfile, datasetName + "_n");
        
        startPosition = 0;
        for (hsize_t i = 0; i < N; ++i)
        {
          HDF5Sequence sequence(NULL, &readIdxArray, &readNameArray, i);
          CuAssertTrue(testCase,
                       sequence.getName() == genName(i, length));
          CuAssertTrue(testCase, 
                       sequence.getStartPosition() == startPosition);
          CuAssertTrue(testCase,
                       sequence.getSequenceLength() == i * 2);
          CuAssertTrue(testCase,
                       sequence.getNumTopSegments() == i * 3);
          CuAssertTrue(testCase,
                       sequence.getNumBottomSegments() == i * 4);
          startPosition += i * 2;
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
