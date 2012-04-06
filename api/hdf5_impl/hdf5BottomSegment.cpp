/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include "hdf5BottomSegment.h"

using namespace std;
using namespace H5;
using namespace hal;

const hsize_t HDF5BottomSegment::genomeIndexOffset = 0;
const hsize_t HDF5BottomSegment::lengthOffset = sizeof(hal_index_t);
const hsize_t HDF5BottomSegment::topIndexOffset = 2 * sizeof(hal_index_t);
const hsize_t HDF5BottomSegment::topOffsetOffset = 3 * sizeof(hal_index_t);
const hsize_t HDF5BottomSegment::parIndexOffset = 3 * sizeof(hal_index_t) + sizeof(hal_offset_t);
const hsize_t HDF5BottomSegment::firstChildOffset = 4 * sizeof(hal_index_t) + sizeof(hal_offset_t);

static const PredType h5index_type = PredType::NATIVE_INT64;
static const PredType h5bool_type = PredType::NATIVE_HBOOL;

H5::CompType HDF5BottomSegment::dataType(hal_size_t numChildren)
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(h5index_type.getSize() == sizeof(hal_index_t));
  assert(h5bool_type.getSize() == sizeof(hal_bool_t));

  H5::CompType dataType;
  dataType.insertMember("genomeIdx", genomeIndexOffset, h5index_type);
  dataType.insertMember("length", lengthOffset, h5index_type);
  dataType.insertMember("topIdx", topIndexOffset, h5index_type);
  dataType.insertMember("topOffset", topOffsetOffset, h5index_type);
  dataType.insertMember("paralogyIdx", parIndexOffset, h5index_type);
  for(hsize_t i = 0; i < numChildren; ++i)
  {
    std::stringstream ss;
    ss << i;
    H5std_string number = ss.str();
    dataType.insertMember("childIdx" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(hal_bool_t)), 
                          h5index_type);
    dataType.insertMember("reverseFlag" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(hal_bool_t)) +
                          sizeof(hal_bool_t), 
                          h5index_type);
  }
  return dataType;
}

