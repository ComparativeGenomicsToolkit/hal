/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/**
 * Test that the basic hdf5 C++ API functionality is working properly
 */

#include <iostream>
#include <string>
#include <H5Cpp.h>
#include "allTests.h"
#include "hdf5Test.h"
extern "C" {
#include "commonC.h"
}

using namespace H5;
using namespace std;

char* fileName = NULL;
int64_t* numbers = NULL;

void hdf5TestTeardown()
{
  delete [] numbers;
  numbers = NULL;
  removeTempFile(fileName);
}

void hdf5TestSetup()
{
  numbers = new int64_t[N];
  for (hsize_t i = 0; i < N; ++i)
  {
    numbers[i] = (int64_t)i;
  }
  fileName = getTempFile();
}

/** Write array of integers (value = index) in h5 file */
void writeNumbers(hsize_t chunkSize)
{
  // write numbers[] to a rawH5ExternalArray file
  H5File file(H5std_string(fileName), H5F_ACC_TRUNC);
  IntType datatype(PredType::NATIVE_HSIZE);
  datatype.setOrder(H5T_ORDER_LE);
  DSetCreatPropList cparms;
  if (chunkSize > N)
     chunkSize = 0;
  if (chunkSize > 0)
  {
    cparms.setDeflate(3);
    cparms.setChunk(1, &chunkSize);
  }
  DataSpace dataspace(1, &N);
  DataSet dataset = file.createDataSet(datasetName, datatype, dataspace, 
                                       cparms);
  if (chunkSize == 0)
  {
    dataset.write(numbers, PredType::NATIVE_HSIZE);
  }
  else
  {
    hsize_t numChunks = N / chunkSize;
    hsize_t lastChunkSize = N % chunkSize != 0 ?  N % chunkSize : chunkSize;
    if (lastChunkSize != chunkSize)
    {
      numChunks++;
    } 
    for (hsize_t i = 0; i < numChunks; ++i)
    {
      hsize_t bufSize = i != numChunks - 1 ? chunkSize : lastChunkSize;
      hsize_t bufStart = i * chunkSize;  
      dataspace.selectHyperslab(H5S_SELECT_SET, &bufSize, &bufStart);
      DataSpace chunkSpace(1, &bufSize);
      dataset.write(numbers + bufStart, PredType::NATIVE_HSIZE, 
                    chunkSpace, dataspace);
    }
  }
  file.close();
}

/** Read an array of integers and verify their values = index */
void checkNumbers(CuTest *testCase)
{
  // read numbers back into memory
  H5File rfile(H5std_string(fileName), H5F_ACC_RDONLY);
  DataSet rdataset = rfile.openDataSet(datasetName);
  H5T_class_t type_class = rdataset.getTypeClass();
  CuAssertTrue(testCase, type_class == H5T_INTEGER);
  IntType rintype = rdataset.getIntType();
    
  DataSpace rdataspace = rdataset.getSpace();
  int rank = rdataspace.getSimpleExtentNdims();
  CuAssertTrue(testCase, rank == 1);

  hsize_t dims_out;
  rdataspace.getSimpleExtentDims(&dims_out, NULL);
  CuAssertTrue(testCase, dims_out == N);

  int64_t* rnumbers = new int64_t[N];
  DataType datatype = rdataset.getDataType();

  DSetCreatPropList cparms = rdataset.getCreatePlist();
  if (cparms.getLayout() == H5D_CHUNKED)
  {
    hsize_t chunkSize;
    cparms.getChunk(1, &chunkSize);
    hsize_t numChunks = N / chunkSize;
    hsize_t lastChunkSize = N % chunkSize != 0 ?  N % chunkSize : chunkSize;
    if (lastChunkSize != chunkSize)
    {
      numChunks++;
    } 
    for (hsize_t i = 0; i < numChunks; ++i)
    {
      hsize_t bufSize = i != numChunks - 1 ? chunkSize : lastChunkSize;
      hsize_t bufStart = i * chunkSize;  
      rdataspace.selectHyperslab(H5S_SELECT_SET, &bufSize, &bufStart);
      DataSpace chunkSpace(1, &bufSize);
      assert(bufStart + bufSize <= N);
      rdataset.read(rnumbers + bufStart, PredType::NATIVE_HSIZE, 
                    chunkSpace, rdataspace);
    }
  }
  else
  {
    rdataset.read(rnumbers, datatype, DataSpace(1, &N), rdataspace);
  }
  for (hsize_t i = 0; i < N; ++i)
  {
    if (rnumbers[i] != numbers[i])
    {
      cerr << i << ": " << rnumbers[i] << " (should be " << numbers[i] << ")"
           << endl;
    }
    CuAssertTrue(testCase, rnumbers[i] == numbers[i]);
  }
  delete [] rnumbers;
}

void hdf5TestBasicFileIO(CuTest *testCase) 
{
  hdf5TestSetup();
  try
  {
    writeNumbers(0);
    checkNumbers(testCase);
  }
  catch(...)
  {
    CuAssertTrue(testCase, 0);
  }
  hdf5TestTeardown();

  hdf5TestSetup();
  try
  {
    writeNumbers(10000);
    checkNumbers(testCase);
  }
  catch(...)
  {
    CuAssertTrue(testCase, 0);
  }
  hdf5TestTeardown();
}

CuSuite* hdf5TestSuite(void) 
{
  CuSuite* suite = CuSuiteNew();
  SUITE_ADD_TEST(suite, hdf5TestBasicFileIO);
  return suite;
}
