/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "hdf5TopSegment.h"
#include "hdf5BottomSegment.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5TopSegment::genomeIndexOffset = 0;
const size_t HDF5TopSegment::lengthOffset = sizeof(hal_index_t);
const size_t HDF5TopSegment::bottomIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5TopSegment::bottomOffsetOffset = bottomIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parIndexOffset = bottomOffsetOffset + sizeof(hal_offset_t);
const size_t HDF5TopSegment::parReversedOffset = parIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parentIndexOffset = parReversedOffset + sizeof(hal_bool_t);
const size_t HDF5TopSegment::parentReversedOffset = parentIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::totalSize = parentReversedOffset + sizeof(hal_bool_t);

HDF5TopSegment::HDF5TopSegment(HDF5Genome* genome,
                               HDF5ExternalArray* array,
                               hal_index_t index) :
  _array(array),
  _index(index),
  _genome(genome)
{
  assert(_index >= 0);
}

HDF5TopSegment::~HDF5TopSegment()
{
  
}

bool HDF5TopSegment::isGapInsertion() const
{
  if (getParentIndex() != NULL_INDEX || 
      getNextParalogyIndex() != NULL_INDEX)
  {
    return false;
  }
  
  const Sequence* sequence = getSequence();
  HDF5Genome* parentGenome =  reinterpret_cast<HDF5Genome*>(
    _genome->getParent());
  hal_index_t seqStartIndex =  sequence->getTopSegmentArrayIndex();
  hal_index_t seqLastIndex = seqStartIndex + 
     (hal_index_t)sequence->getNumTopSegments() - 1;

  hal_index_t leftParentIndex = NULL_INDEX;
  const Sequence* leftParentSequence = NULL;
  hal_index_t rightParentIndex = NULL_INDEX;
  const Sequence* rightParentSequence = NULL;

  // walk up to parent of left neighbour
  if (_index != seqStartIndex)
  {
    HDF5TopSegment leftSeg(_genome, _array, _index - 1);
    if (leftSeg.getParentIndex() != NULL_INDEX)
    {
      HDF5BottomSegment leftPar(parentGenome, &parentGenome->_bottomArray,
                                leftSeg.getParentIndex());
      leftParentIndex = leftPar.getArrayIndex();
      leftParentSequence = leftPar.getSequence();
    }
  }

  // wlak up to parent of right neighbour
  if (_index != seqLastIndex)
  {
    HDF5TopSegment rightSeg(_genome, _array, _index + 1);
    if (rightSeg.getParentIndex() != NULL_INDEX)
    {
      HDF5BottomSegment rightPar(parentGenome, &parentGenome->_bottomArray,
                                rightSeg.getParentIndex());
      rightParentIndex = rightPar.getArrayIndex();
      rightParentSequence = rightPar.getSequence();
    }
  }

  // case 1) gap insertion inside a sequence
  if (leftParentIndex != NULL_INDEX && rightParentIndex != NULL_INDEX &&
      abs(rightParentIndex - leftParentIndex) == 1 &&
      leftParentSequence == rightParentSequence)
  {
    return true;
  }

  // case 2) gap insertion at beginning of sequence
  if (leftParentIndex == NULL_INDEX && rightParentIndex != NULL_INDEX &&
      rightParentIndex == rightParentSequence->getBottomSegmentArrayIndex())
  {
    return true;
  }
  
  // case 3) gap insertion at end of sequence
  if (rightParentIndex == NULL_INDEX && leftParentIndex != NULL_INDEX &&
      leftParentIndex == leftParentSequence->getBottomSegmentArrayIndex() +
      (hal_index_t)leftParentSequence->getNumBottomSegments() - 1)
  {
    return true;
  }

  return false;
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
  dataType.insertMember("paralogyReversed", parReversedOffset, PredType::NATIVE_CHAR);
  dataType.insertMember("parentIdx", parentIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("reverseFlag", parentReversedOffset, PredType::NATIVE_CHAR);

  return dataType;
}

