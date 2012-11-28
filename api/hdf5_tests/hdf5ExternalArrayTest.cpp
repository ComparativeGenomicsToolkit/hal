/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * Test the external HDF5 array wrapper
 */

#include <iostream>
#include <string>
#include <H5Cpp.h>
#include "allTests.h"
#include "hdf5ExternalArray.h"
#include "hdf5Test.h"
extern "C" {
#include "commonC.h"
}

using namespace H5;
using namespace hal;
using namespace std;

static const hsize_t chunkSizes[] = {0, 4, N/10, N/5, N/2, N, 2*N};
static const hsize_t numSizes = 7;

static void teardown()
{
  hdf5TestTeardown();
}

static void setup()
{
  hdf5TestSetup();
}

void hdf5ExternalArrayTestCreate(CuTest *testCase)
{
  for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
  {
    hsize_t chunkSize = chunkSizes[chunkIdx];
    setup();
    try 
    {
      IntType datatype(PredType::NATIVE_HSIZE);

      H5File file(H5std_string(fileName), H5F_ACC_TRUNC);

      HDF5ExternalArray myArray;
      DSetCreatPropList cparms;
      if (chunkSize > 0)
      {
        cparms.setDeflate(2);
        cparms.setChunk(1, &chunkSize);
      }
      myArray.create(&file, datasetName, datatype, N, &cparms);
      for (hsize_t i = 0; i < N; ++i)
      {
        hsize_t* block = reinterpret_cast<hsize_t*>(myArray.getUpdate(i));
        *block = i;
      }
      myArray.write();
      file.flush(H5F_SCOPE_LOCAL);
      file.close();
      checkNumbers(testCase);
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

void hdf5ExternalArrayTestLoad(CuTest *testCase)
{
  for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
  {
    hsize_t chunkSize = chunkSizes[chunkIdx];
    setup();
    try 
    {
      writeNumbers(chunkSize);
      
      H5File file(H5std_string(fileName), H5F_ACC_RDONLY);
      HDF5ExternalArray myArray;
      myArray.load(&file, datasetName);

      for (hsize_t i = 0; i < N; ++i)
      {
        const int64_t* val = reinterpret_cast<const int64_t*>(myArray.get(i));
        CuAssertTrue(testCase, *val == numbers[i]);
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
  }
}

void hdf5ExternalArrayTestCompression(CuTest *testCase)
{
  for (hsize_t chunkIdx = 0; chunkIdx < numSizes; ++chunkIdx)
  {
    hsize_t chunkSize = chunkSizes[chunkIdx];
    setup();
    try 
    {
      IntType datatype(PredType::NATIVE_HSIZE);
      H5File file(H5std_string(fileName), H5F_ACC_TRUNC);
      HDF5ExternalArray myArray;
      DSetCreatPropList cparms;
      if (chunkSize > 0)
      {
        cparms.setDeflate(9);
        cparms.setChunk(1, &chunkSize);
      }
      myArray.create(&file, datasetName, datatype, N, &cparms);
      for (hsize_t i = 0; i < N; ++i)
      {
        hsize_t* block = reinterpret_cast<hsize_t*>(myArray.getUpdate(i));
        *block = i;
      }
      myArray.write();
      file.flush(H5F_SCOPE_LOCAL);
      file.close();

      H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
      HDF5ExternalArray myrArray;
      myrArray.load(&rfile, datasetName);
      
      for (hsize_t i = 0; i < N; ++i)
      {
        const int64_t* val = reinterpret_cast<const int64_t*>(myrArray.get(i));
        CuAssertTrue(testCase, *val == numbers[i]);
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
  }
}

CuSuite* hdf5ExternalArrayTestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, hdf5ExternalArrayTestCreate);
  SUITE_ADD_TEST(suite, hdf5ExternalArrayTestLoad);
  SUITE_ADD_TEST(suite, hdf5ExternalArrayTestCompression);
  return suite;
}
