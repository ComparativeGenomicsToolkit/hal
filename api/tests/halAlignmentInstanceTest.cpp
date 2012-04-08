/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halAlignmentInstanceTest.h"
#include "halAlignmentInstance.h"

using namespace std;
using namespace hal;

vector<AlignmentPtr> getTestAlignmentInstances()
{
  vector<AlignmentPtr> testInstances;

  // DEFAULT HDF5
  testInstances.push_back(hdf5AlignmentInstance());
  
  // TODO : CHUNKING, CACHING, COMPRESSION

  return testInstances;
}
