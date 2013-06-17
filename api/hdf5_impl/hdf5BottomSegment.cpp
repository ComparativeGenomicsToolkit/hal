/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <cstdlib>
#include "hdf5BottomSegment.h"
#include "hdf5TopSegment.h"
#include "hdf5DNAIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5BottomSegment::genomeIndexOffset = 0;
const size_t HDF5BottomSegment::lengthOffset = sizeof(hal_index_t);
const size_t HDF5BottomSegment::topIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5BottomSegment::firstChildOffset = topIndexOffset + sizeof(hal_index_t);
const size_t HDF5BottomSegment::totalSize(hal_size_t numChildren)
{
  return firstChildOffset + numChildren * (sizeof(hal_index_t) + sizeof(bool));
}

HDF5BottomSegment::HDF5BottomSegment(HDF5Genome* genome,
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

hal_size_t HDF5BottomSegment::numChildrenFromDataType(
  const H5::DataType& dataType)
{
  return (dataType.getSize() - firstChildOffset) / 
     (sizeof(hal_index_t) + sizeof(bool));
}

hal_offset_t HDF5BottomSegment::getTopParseOffset() const
{
  assert(_index >= 0);
  hal_offset_t offset = 0;
  hal_index_t topIndex = getTopParseIndex();
  if (topIndex != NULL_INDEX)
  {
    HDF5TopSegment ts(_genome, &_genome->_topArray, topIndex);
    assert(ts.getStartPosition() <= getStartPosition());
    assert((hal_index_t)(ts.getStartPosition() + ts.getLength()) 
           >= getStartPosition());
    offset = getStartPosition() - ts.getStartPosition();
  }
  return offset;
}

void HDF5BottomSegment::setCoordinates(hal_index_t startPos, hal_size_t length)
{
  assert(_index >= 0);
  if (_genome && (startPos >= (hal_index_t)_genome->_totalSequenceLength || 
                  startPos + length > _genome->_totalSequenceLength))
  {
    throw hal_exception("Trying to set bottom segment coordinate out of range");
  }
  
  _array->setValue((hsize_t)_index, genomeIndexOffset, startPos);
  _array->setValue(_index + 1, genomeIndexOffset, startPos + length);
}

void HDF5BottomSegment::getString(string& outString) const
{
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_genome), getStartPosition());
  di.readString(outString, getLength()); 
}

bool HDF5BottomSegment::isMissingData(double nThreshold) const
{
  if (nThreshold >= 1.0)
  {
    return false;
  }  
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_genome), getStartPosition());
  size_t length = getLength();
  size_t maxNs = nThreshold * (double)length;
  size_t Ns = 0;
  char c;
  for (size_t i = 0; i < length; ++i, di.toRight())
  {
    c = di.getChar();
    if (c == 'N' || c == 'n')
    {
      ++Ns;
    }
    if (Ns > maxNs)
    {
      return true;
    }
    if ((length - i) < (maxNs - Ns))
    {
      break;
    }
  }
  return false;
}

void HDF5BottomSegment::print(std::ostream& os) const
{
  os << "HDF5 Bottom Segment";
}

// HDF5 SPECIFIC
H5::CompType HDF5BottomSegment::dataType(hal_size_t numChildren)
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
  assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(bool));

  H5::CompType dataType(totalSize(numChildren));
  dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("length", lengthOffset, PredType::NATIVE_HSIZE);
  dataType.insertMember("topIdx", topIndexOffset, PredType::NATIVE_INT64);
  for(hsize_t i = 0; i < numChildren; ++i)
  {
    std::stringstream ss;
    ss << i;
    H5std_string number = ss.str();
    dataType.insertMember("childIdx" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(bool)), 
                          PredType::NATIVE_INT64);
    dataType.insertMember("reverseFlag" + number, firstChildOffset + 
                          i * (sizeof(hal_index_t) + sizeof(bool)) +
                          sizeof(hal_index_t), 
                          PredType::NATIVE_CHAR);
  }
  return dataType;
}

