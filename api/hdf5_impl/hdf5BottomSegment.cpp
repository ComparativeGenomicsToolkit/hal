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

const size_t HDF5BottomSegment::genomeIndexOffset = 0;
const size_t HDF5BottomSegment::lengthOffset = sizeof(hal_index_t);
const size_t HDF5BottomSegment::topIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5BottomSegment::topOffsetOffset = topIndexOffset + sizeof(hal_index_t);
const size_t HDF5BottomSegment::parIndexOffset = topOffsetOffset + sizeof(hal_offset_t);
const size_t HDF5BottomSegment::firstChildOffset = parIndexOffset + sizeof(hal_index_t);
const size_t HDF5BottomSegment::totalSize(hal_size_t numChildren)
{
  return firstChildOffset + numChildren * (sizeof(hal_index_t) + sizeof(hal_bool_t));
}

HDF5BottomSegment::HDF5BottomSegment(GenomePtr genome,
                                     HDF5ExternalArray* array,
                                     hal_index_t index) :
  _array(array),
  _index(index),
  _genome(genome)
{

}

HDF5BottomSegment::~HDF5BottomSegment()
{
  
}

H5::CompType HDF5BottomSegment::dataType(hal_size_t numChildren)
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
  assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(hal_bool_t));

  H5::CompType dataType(totalSize(numChildren));
  dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("length", lengthOffset, PredType::NATIVE_HSIZE);
  dataType.insertMember("topIdx", topIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("topOffset", topOffsetOffset, PredType::NATIVE_UINT64);
  dataType.insertMember("paralogyIdx", parIndexOffset, PredType::NATIVE_INT64);
  for(hsize_t i = 0; i < numChildren; ++i)
  {
    std::stringstream ss;
    ss << i;
    H5std_string number = ss.str();
    dataType.insertMember("childIdx" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(hal_bool_t)), 
                          PredType::NATIVE_INT64);
    dataType.insertMember("reverseFlag" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(hal_bool_t)) +
                          sizeof(hal_index_t), 
                          PredType::NATIVE_CHAR);
  }
  return dataType;
}

