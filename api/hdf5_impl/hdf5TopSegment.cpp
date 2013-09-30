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
#include "hdf5DNAIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5TopSegment::genomeIndexOffset = 0;
const size_t HDF5TopSegment::bottomIndexOffset = sizeof(hal_index_t);
const size_t HDF5TopSegment::parIndexOffset = bottomIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parentIndexOffset = parIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::parentReversedOffset = parentIndexOffset + sizeof(hal_index_t);
const size_t HDF5TopSegment::totalSize = parentReversedOffset + sizeof(bool);

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

void HDF5TopSegment::setCoordinates(hal_index_t startPos, hal_size_t length)
{
  if (_genome && (startPos >= (hal_index_t)_genome->_totalSequenceLength || 
                  startPos + length > _genome->_totalSequenceLength))
  {
    throw hal_exception("Trying to set top segment coordinate out of range");
  }
      
  _array->setValue(_index, genomeIndexOffset, startPos);
  _array->setValue(_index + 1, genomeIndexOffset, startPos + length);
}
   
hal_offset_t HDF5TopSegment::getBottomParseOffset() const
{
  assert(_index >= 0);
  hal_offset_t offset = 0;
  hal_index_t bottomIndex = getBottomParseIndex();
  if (bottomIndex != NULL_INDEX)
  {
    HDF5BottomSegment bs(_genome, &_genome->_bottomArray, bottomIndex);
    assert(bs.getStartPosition() <= getStartPosition());
    assert((hal_index_t)(bs.getStartPosition() + bs.getLength()) 
           >= getStartPosition());
    offset = getStartPosition() - bs.getStartPosition();
  }
  return offset;
}

void HDF5TopSegment::getString(std::string& outString) const
{
  HDF5DNAIterator di(const_cast<HDF5Genome*>(_genome), getStartPosition());
  di.readString(outString, getLength()); 
}

bool HDF5TopSegment::isMissingData(double nThreshold) const
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

bool HDF5TopSegment::isCanonicalParalog() const
{
  bool isCanon = false;
  if (hasParent())
  {
    HDF5Genome* parGenome = 
       const_cast <HDF5Genome*>(
         dynamic_cast<const HDF5Genome*>(_genome->getParent()));

    HDF5BottomSegment parent(parGenome, 
                             &parGenome->_bottomArray,
                             getParentIndex());
    hal_index_t childGenomeIndex = parGenome->getChildIndex(_genome);
    isCanon = parent.getChildIndex(childGenomeIndex) == _index;
  }
  return isCanon;
}

void HDF5TopSegment::print(std::ostream& os) const
{
  os << "HDF5 Top Segment";
}

// HDF5 SPECIFIC
H5::CompType HDF5TopSegment::dataType()
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
  assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(bool));

  H5::CompType dataType(totalSize);
  dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("bottomIdx", bottomIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("paralogyIdx", parIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("parentIdx", parentIndexOffset, PredType::NATIVE_INT64);
  dataType.insertMember("reverseFlag", parentReversedOffset, PredType::NATIVE_CHAR);

  return dataType;
}
