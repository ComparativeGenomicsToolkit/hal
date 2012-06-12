/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hdf5DNAIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

HDF5DNAIterator::HDF5DNAIterator(HDF5Genome* genome, hal_index_t index) :
  _index(index),
  _genome(genome),
  _reversed(false)
{

}

HDF5DNAIterator::~HDF5DNAIterator()
{

}
