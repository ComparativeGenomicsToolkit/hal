/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include "hdf5TopSegment.h"

using namespace std;
using namespace H5;
using namespace hal;

const hsize_t HDF5TopSegment::genomeIndexOffset = 0;
const hsize_t HDF5TopSegment::lengthOffset = sizeof(hal_index_t);
const hsize_t HDF5TopSegment::bottomIndexOffset = 2 * sizeof(hal_index_t);
const hsize_t HDF5TopSegment::bottomOffsetOffset = 3 * sizeof(hal_index_t);
const hsize_t HDF5TopSegment::parIndexOffset = 3 * sizeof(hal_index_t) + sizeof(hal_offset_t);
const hsize_t HDF5TopSegment::parentIndexOffset = 4 * sizeof(hal_index_t) + sizeof(hal_offset_t);
const hsize_t HDF5TopSegment::parentReversedOffset = 5 * sizeof(hal_index_t) + sizeof(hal_offset_t);

static const PredType h5index_type = PredType::NATIVE_INT64;
static const PredType h5bool_type = PredType::NATIVE_HBOOL;

H5::CompType HDF5TopSegment::dataType()
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(h5index_type.getSize() == sizeof(hal_index_t));
  assert(h5bool_type.getSize() == sizeof(hal_bool_t));

  H5::CompType dataType;
  dataType.insertMember("genomeIdx", genomeIndexOffset, h5index_type);
  dataType.insertMember("length", lengthOffset, h5index_type);
  dataType.insertMember("bottomIdx", bottomIndexOffset, h5index_type);
  dataType.insertMember("bottomOffset", bottomOffsetOffset, h5index_type);
  dataType.insertMember("paralogyIdx", parIndexOffset, h5index_type);
  dataType.insertMember("parentIdx", parentIndexOffset, h5index_type);
  dataType.insertMember("reverseFlag", parentReversedOffset, h5bool_type);

  return dataType;
}

