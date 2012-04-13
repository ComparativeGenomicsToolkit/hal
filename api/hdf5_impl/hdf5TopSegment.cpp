/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include "hdf5TopSegment.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5TopSegment::genomeIndexOffset = 0;
const size_t HDF5TopSegment::lengthOffset = sizeof(hal_index_t);
const size_t HDF5TopSegment::bottomIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5TopSegment::bottomOffsetOffset = bottomIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parIndexOffset = bottomOffsetOffset + sizeof(hal_offset_t);
const size_t HDF5TopSegment::parentIndexOffset = parIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parentReversedOffset = parentIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::totalSize = parentReversedOffset + sizeof(hal_bool_t);

HDF5TopSegment::HDF5TopSegment(HDF5Genome* genome,
                               HDF5ExternalArray* array,
                               hal_index_t index) :
  _array(array),
  _index(index),
  _genome(genome)
{
  
}

HDF5TopSegment::~HDF5TopSegment()
{
  
}

H5::CompType HDF5TopSegment::dataType()
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
  assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(hal_bool_t));

  H5::CompType dataType(totalSize);
  dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("length", lengthOffset, PredType::NATIVE_HSIZE); 
  dataType.insertMember("bottomIdx", bottomIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("bottomOffset", bottomOffsetOffset, PredType::NATIVE_UINT64);
  dataType.insertMember("paralogyIdx", parIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("parentIdx", parentIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("reverseFlag", parentReversedOffset, PredType::NATIVE_CHAR);

  return dataType;
}

