/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5TEST_H
#define _HDF5TEST_H

/** some functionality shared by tests that rely on basic
 * hdf5 stuff
 */

static const hsize_t N = 500000;
static const std::string datasetName("name");
extern char* fileName;
extern int64_t* numbers;

void hdf5TestTeardown();
void hdf5TestSetup();
void writeNumbers(hsize_t chunkSize);
void checkNumbers(CuTest *testCase);

#endif
