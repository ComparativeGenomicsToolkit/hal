/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include "hdf5Sequence.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5Sequence::startOffset = 0;
const size_t HDF5Sequence::lengthOffset = sizeof(hal_size_t);
const size_t HDF5Sequence::numTopSegmentsOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::numBottomSegmentsOffset = numTopSegmentsOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::nameOffset = numBottomSegmentsOffset + sizeof(hal_size_t);

HDF5Sequence::HDF5Sequence(HDF5Genome* genome,
                           HDF5ExternalArray* array,
                           hal_index_t index) :
  _array(array),
  _index(index),
  _genome(genome)
{

}

HDF5Sequence::~HDF5Sequence()
{
  
}

H5::CompType HDF5Sequence::dataType(hal_size_t maxNameLength)
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(char));
  StrType strType(PredType::NATIVE_CHAR, (maxNameLength + 1) * sizeof(char));
  CompType dataType(nameOffset + strType.getSize());
  dataType.insertMember("start", startOffset, PredType::NATIVE_HSIZE);
  dataType.insertMember("length", lengthOffset, PredType::NATIVE_HSIZE);
  dataType.insertMember("numSequences", numTopSegmentsOffset, 
                        PredType::NATIVE_HSIZE);
  dataType.insertMember("numBottomSegments", numBottomSegmentsOffset, 
                        PredType::NATIVE_HSIZE);
  dataType.insertMember("name", nameOffset, strType);
  return dataType;
}

// SEQUENCE INTERFACE
string HDF5Sequence::getName() const
{
  return _array->get(_index) + nameOffset;
}

const Genome* HDF5Sequence::getGenome() const
{
  return _genome;
}

Genome* HDF5Sequence::getGenome()
{
  return _genome;
}

hal_size_t HDF5Sequence::getStartPosition() const
{
  return _array->getValue<hal_size_t>(_index, startOffset);
}

hal_index_t HDF5Sequence::getArrayIndex() const
{
  return _index;
}

// SEGMENTED SEQUENCE INTERFACE

hal_size_t HDF5Sequence::getSequenceLength() const
{
  return _array->getValue<hal_size_t>(_index, lengthOffset);
}

hal_size_t HDF5Sequence::getNumTopSegments() const
{
  return _array->getValue<hal_size_t>(_index, numTopSegmentsOffset);
}

hal_size_t HDF5Sequence::getNumBottomSegments() const
{
  return _array->getValue<hal_size_t>(_index, numBottomSegmentsOffset);
}

TopSegmentIteratorPtr HDF5Sequence::getTopSegmentIterator(
  hal_index_t position)
{
  
  return TopSegmentIteratorPtr();
}

TopSegmentIteratorConstPtr HDF5Sequence::getTopSegmentIterator(
  hal_index_t position) const
{
  return TopSegmentIteratorConstPtr();
}

BottomSegmentIteratorPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position)
{
  return BottomSegmentIteratorPtr();
}

BottomSegmentIteratorConstPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position) const
{
  return BottomSegmentIteratorConstPtr();
}

void HDF5Sequence::getString(std::string& outString) const
{

}

void HDF5Sequence::setString(const std::string& inString)
{

}

void HDF5Sequence::getSubString(std::string& outString, hal_size_t start,
                                hal_size_t length) const
{

}

void HDF5Sequence::setSubString(const std::string& intString, 
                                hal_size_t start,
                                hal_size_t length)
{

}

// LOCAL

void HDF5Sequence::set(hal_size_t startPosition, 
                       const Sequence::Info& sequenceInfo)
{
  _array->setValue(_index, startOffset, startPosition);
  _array->setValue(_index, lengthOffset, sequenceInfo._length);
  _array->setValue(_index, numTopSegmentsOffset, sequenceInfo._numTopSegments);
  _array->setValue(_index, numBottomSegmentsOffset,
                   sequenceInfo._numBottomSegments);
  char* arrayBuffer = _array->getUpdate(_index) + nameOffset;
  strcpy(arrayBuffer, sequenceInfo._name.c_str());
}

