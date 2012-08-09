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

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5BottomSegment::genomeIndexOffset = 0;
const size_t HDF5BottomSegment::lengthOffset = sizeof(hal_index_t);
const size_t HDF5BottomSegment::topIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5BottomSegment::topOffsetOffset = topIndexOffset + sizeof(hal_index_t);
const size_t HDF5BottomSegment::parIndexOffset = topOffsetOffset + sizeof(hal_offset_t);
const size_t HDF5BottomSegment::parReversedOffset = parIndexOffset + sizeof(hal_index_t);
const size_t HDF5BottomSegment::firstChildOffset = parReversedOffset + sizeof(hal_bool_t);
const size_t HDF5BottomSegment::totalSize(hal_size_t numChildren)
{
  return firstChildOffset + numChildren * (sizeof(hal_index_t) + sizeof(hal_bool_t));
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

bool HDF5BottomSegment::isGapDeletion(hal_size_t i) const
{
  if (getChildIndex(i) != NULL_INDEX ||
      getNextParalogyIndex() != NULL_INDEX)
  {
    return false;
  }

  HDF5Genome* childGenome =  reinterpret_cast<HDF5Genome*>(
    _genome->getChild(i));

  hal_index_t leftChildIndex = NULL_INDEX;
  const Sequence* leftChildSequence = NULL;
  hal_index_t rightChildIndex = NULL_INDEX;
  const Sequence* rightChildSequence = NULL;

  // walk down to child of left neighbour
  if (isFirst() == false)
  {
    leftChildIndex = getLeftChildIndex(i);
    if (leftChildIndex != NULL_INDEX)
    {
      HDF5TopSegment leftChild(childGenome, &childGenome->_topArray,
                               leftChildIndex);
      leftChildSequence = leftChild.getSequence();
    }
  }

  // wlak down to child of right neighbour
  if (isLast() == false)
  {
    rightChildIndex = getRightChildIndex(i);
    if (rightChildIndex != NULL_INDEX)
    {
      HDF5TopSegment rightChild(childGenome, &childGenome->_topArray,
                                rightChildIndex);
      rightChildSequence = rightChild.getSequence();
    }
  }

  // case 1) gap deletion inside a sequence
  if (leftChildIndex != NULL_INDEX && rightChildIndex != NULL_INDEX &&
      abs(rightChildIndex - leftChildIndex) == 1 &&
      leftChildSequence == rightChildSequence)
  {
    return true;
  }

  // case 2) gap deletion at beginning of sequence
  if (leftChildIndex == NULL_INDEX && rightChildIndex != NULL_INDEX &&
      rightChildIndex == rightChildSequence->getTopSegmentArrayIndex())
  {
    return true;
  }
  
  // case 3) gap insertion at end of sequence
  if (rightChildIndex == NULL_INDEX && leftChildIndex != NULL_INDEX &&
      leftChildIndex == leftChildSequence->getTopSegmentArrayIndex() +
      (hal_index_t)leftChildSequence->getNumTopSegments() - 1)
  {
    return true;
  }

  return false;
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
  dataType.insertMember("paralogyReverseFlag", parReversedOffset, 
                        PredType::NATIVE_CHAR);
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

hal_size_t HDF5BottomSegment::numChildrenFromDataType(
  const H5::DataType& dataType)
{
  return (dataType.getSize() - firstChildOffset) / 
     (sizeof(hal_index_t) + sizeof(hal_bool_t));
}
