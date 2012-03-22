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
extern "C" {
#include "commonC.h"
}

using namespace H5;

static char* fileName = NULL;
static const hsize_t N = 1000;
static int64_t numbers[N];

static void hdf5TestTeardown()
{
	removeTempFile(fileName);
}

static void hdf5TestSetup()
{
	for (hsize_t i = 0; i < N; ++i)
	{
		numbers[i] = (int64_t)i;
	}
	fileName = getTempFile();
}

void testBasicFileIO(CuTest *testCase) 
{
    hdf5TestSetup();
	 try
	 {
		 // write numbers[] to a hdf5 file
		 H5File file(H5std_string(fileName), H5F_ACC_TRUNC);
		 DataSpace dataspace(1, &N);
		 IntType datatype(PredType::NATIVE_LONG);
		 datatype.setOrder(H5T_ORDER_LE);
		 H5std_string datasetName = "dataset";
		 DataSet dataset = file.createDataSet(datasetName, datatype, dataspace);
		 dataset.write(numbers, PredType::NATIVE_LONG);
		 file.close();

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

		 int64_t rnumbers[N];
		 rdataset.read(rnumbers, datatype);
		 for (hsize_t i = 0; i < N; ++i)
		 {
			 CuAssertTrue(testCase, rnumbers[i] == numbers[i]);
		 }
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
    SUITE_ADD_TEST(suite, testBasicFileIO);
    return suite;
}
