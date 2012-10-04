/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include "hdf5Sequence.h"
#include "hdf5DNAIterator.h"
#include "hdf5TopSegmentIterator.h"
#include "hdf5BottomSegmentIterator.h"
#include "defaultColumnIterator.h"
#include "defaultRearrangement.h"
#include "defaultGappedTopSegmentIterator.h"
#include "defaultGappedBottomSegmentIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5Sequence::startOffset = 0;
const size_t HDF5Sequence::lengthOffset = sizeof(hal_size_t);
const size_t HDF5Sequence::numTopSegmentsOffset = lengthOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::numBottomSegmentsOffset = numTopSegmentsOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::topSegmentArrayIndexOffset = numBottomSegmentsOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::bottomSegmentArrayIndexOffset = topSegmentArrayIndexOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::nameOffset = bottomSegmentArrayIndexOffset + sizeof(hal_size_t);

HDF5Sequence::HDF5Sequence(HDF5Genome* genome,
                           HDF5ExternalArray* array,
                           hal_index_t index) :
  _array(array),
  _index(index),
  _genome(genome),
  _cacheIndex(NULL_INDEX),
  _startCache(NULL_INDEX),
  _lengthCache(0)
{

}

HDF5Sequence::~HDF5Sequence()
{
  
}

// Some tools (ie maf export) access the same sequence
// info over and over again deep inside a loop.  So I hack
// in a little cache, mostly to avoid copying the same string
// (name) out of hdf5.  
inline void HDF5Sequence::refreshCache() const
{
  if (_index != _cacheIndex)
  {
    _nameCache = _array->get(_index) + nameOffset;
    _startCache = _array->getValue<hal_size_t>(_index, startOffset);
    _lengthCache = _array->getValue<hal_size_t>(_index, lengthOffset);
    _cacheIndex = _index;
  }
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
  dataType.insertMember("topSegmentArrayIndexOffset", topSegmentArrayIndexOffset, 
                        PredType::NATIVE_HSIZE);
  dataType.insertMember("bottomSegmentArrayIndexOffset", bottomSegmentArrayIndexOffset, 
                        PredType::NATIVE_HSIZE);
  dataType.insertMember("name", nameOffset, strType);
  return dataType;
}

// SEQUENCE INTERFACE
const string& HDF5Sequence::getName() const
{
  refreshCache();
  return _nameCache;
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
  refreshCache();
  return _startCache;
}

hal_index_t HDF5Sequence::getArrayIndex() const
{
  return _index;
}

hal_index_t HDF5Sequence::getTopSegmentArrayIndex() const
{
  return (hal_index_t)
     _array->getValue<hal_size_t>(_index, topSegmentArrayIndexOffset);
}

hal_index_t HDF5Sequence::getBottomSegmentArrayIndex() const
{
  return (hal_index_t)
     _array->getValue<hal_size_t>(_index, bottomSegmentArrayIndexOffset);
}

// SEGMENTED SEQUENCE INTERFACE

hal_size_t HDF5Sequence::getSequenceLength() const
{
  refreshCache();
  return _lengthCache;
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
  hal_size_t idx = position + getTopSegmentArrayIndex();
  HDF5TopSegmentIterator* newIt = new HDF5TopSegmentIterator(_genome, idx);
  return TopSegmentIteratorPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Sequence::getTopSegmentIterator(
  hal_index_t position) const
{
  hal_size_t idx = position + getTopSegmentArrayIndex();
  HDF5Genome* genome = const_cast<HDF5Genome*>(_genome);
  const HDF5TopSegmentIterator* newIt = 
     new HDF5TopSegmentIterator(genome, idx);
  return TopSegmentIteratorConstPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Sequence::getTopSegmentEndIterator() const
{
  return getTopSegmentIterator(getTopSegmentArrayIndex() + getNumTopSegments());
}

BottomSegmentIteratorPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position)
{
  hal_size_t idx = position + getBottomSegmentArrayIndex();
  HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(_genome, idx);
  return BottomSegmentIteratorPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position) const
{
  hal_size_t idx = position + getBottomSegmentArrayIndex();
  HDF5Genome* genome = const_cast<HDF5Genome*>(_genome);
  const HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(genome, idx);
  return BottomSegmentIteratorConstPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Sequence::getBottomSegmentEndIterator() const
{
  return getBottomSegmentIterator(getBottomSegmentArrayIndex() +
                                  getNumBottomSegments());
}

DNAIteratorPtr HDF5Sequence::getDNAIterator(hal_index_t position)
{
  hal_size_t idx = position + getStartPosition();
  HDF5DNAIterator* newIt = new HDF5DNAIterator(_genome, idx);
  return DNAIteratorPtr(newIt);
}

DNAIteratorConstPtr HDF5Sequence::getDNAIterator(hal_index_t position) const
{
  hal_size_t idx = position + getStartPosition();
  HDF5Genome* genome = const_cast<HDF5Genome*>(_genome);
  const HDF5DNAIterator* newIt = new HDF5DNAIterator(genome, idx);
  return DNAIteratorConstPtr(newIt);
}

DNAIteratorConstPtr HDF5Sequence::getDNAEndIterator() const
{
  return getDNAIterator(getStartPosition() + getSequenceLength());
}

ColumnIteratorConstPtr HDF5Sequence::getColumnIterator(
  const Genome* root, hal_size_t maxInsertLength, hal_index_t position,
  hal_index_t lastPosition, bool noDupes) const
{
  hal_index_t idx = (hal_index_t)(position + getStartPosition());
  hal_index_t lastIdx;
  if (lastPosition == NULL_INDEX)
  {
    lastIdx = (hal_index_t)(getStartPosition() + getSequenceLength() - 1);
  }
  else
  {
    lastIdx = (hal_index_t)(lastPosition + getStartPosition());
  }
  if (position < 0 || 
      lastPosition >= (hal_index_t)(getStartPosition() + getSequenceLength()))
  {
    stringstream ss;
    ss << "HDF5Sequence::getColumnIteratorsetString: input indices "
       << "(" << position << ", " << lastPosition << ") out of bounds";
    throw hal_exception(ss.str());
  }
  const DefaultColumnIterator* newIt = 
     new DefaultColumnIterator(getGenome(), root, idx, lastIdx, maxInsertLength, 
                               noDupes);
  return ColumnIteratorConstPtr(newIt);
}

void HDF5Sequence::getString(std::string& outString) const
{
  getSubString(outString, 0, getSequenceLength());
}

void HDF5Sequence::setString(const std::string& inString)
{
  setSubString(inString, 0, getSequenceLength());
}

void HDF5Sequence::getSubString(std::string& outString, hal_size_t start,
                                hal_size_t length) const
{
  hal_size_t idx = start + getStartPosition();
  outString.resize(length);
  HDF5DNAIterator dnaIt(const_cast<HDF5Genome*>(_genome), idx);
  dnaIt.readString(outString, length);
}

void HDF5Sequence::setSubString(const std::string& inString, 
                                hal_size_t start,
                                hal_size_t length)
{
  if (length != inString.length())
  {
    stringstream ss;
    ss << "setString: input string of length " << inString.length()
       << " has length different from target string in sequence " << getName() 
       << " which is of length " << length;
    throw hal_exception(ss.str());
  }
  hal_size_t idx = start + getStartPosition();
  HDF5DNAIterator dnaIt(_genome, idx);
  dnaIt.writeString(inString, length);
}

RearrangementPtr HDF5Sequence::getRearrangement(hal_index_t position) const
{
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(position);  
  DefaultRearrangement* rea = new DefaultRearrangement(getGenome());
  rea->identifyFromLeftBreakpoint(top);
  return RearrangementPtr(rea);
}

GappedTopSegmentIteratorConstPtr HDF5Sequence::getGappedTopSegmentIterator(
  hal_index_t i, hal_size_t gapThreshold, bool atomic) const
{
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(i);  
  DefaultGappedTopSegmentIterator* gt = 
     new DefaultGappedTopSegmentIterator(top, gapThreshold, atomic);
  return GappedTopSegmentIteratorConstPtr(gt);
}

GappedBottomSegmentIteratorConstPtr 
HDF5Sequence::getGappedBottomSegmentIterator(
  hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
  bool atomic) const
{
  BottomSegmentIteratorConstPtr bot = getBottomSegmentIterator(i);  
  DefaultGappedBottomSegmentIterator* gb = 
     new DefaultGappedBottomSegmentIterator(bot, childIdx, gapThreshold, 
                                            atomic);
  return GappedBottomSegmentIteratorConstPtr(gb);
}
// LOCAL

void HDF5Sequence::set(hal_size_t startPosition, 
                       const Sequence::Info& sequenceInfo,
                       hal_size_t topSegmentStartIndex,
                       hal_size_t bottomSegmentStartIndex)
{
  _array->setValue(_index, startOffset, startPosition);
  _array->setValue(_index, lengthOffset, sequenceInfo._length);
  _array->setValue(_index, numTopSegmentsOffset, sequenceInfo._numTopSegments);
  _array->setValue(_index, numBottomSegmentsOffset,
                   sequenceInfo._numBottomSegments);
  _array->setValue(_index, topSegmentArrayIndexOffset, topSegmentStartIndex);
  _array->setValue(_index, bottomSegmentArrayIndexOffset, 
                   bottomSegmentStartIndex);
  char* arrayBuffer = _array->getUpdate(_index) + nameOffset;
  strcpy(arrayBuffer, sequenceInfo._name.c_str());
  // keep cache for frequently used queries. 
  _cacheIndex = NULL_INDEX;
  refreshCache();
}

void HDF5Sequence::setNumTopSegments(hal_size_t numTopSegments)
{
  _array->setValue(_index, numTopSegmentsOffset, numTopSegments);
}

void HDF5Sequence::setNumBottomSegments(hal_size_t numBottomSegments)
{
  _array->setValue(_index, numBottomSegmentsOffset, numBottomSegments);
}

void HDF5Sequence::setTopSegmentArrayIndex(hal_size_t topIndex)
{
  _array->setValue(_index, topSegmentArrayIndexOffset, topIndex);
}

void HDF5Sequence::setBottomSegmentArrayIndex(hal_size_t bottomIndex)
{
  _array->setValue(_index, bottomSegmentArrayIndexOffset, bottomIndex);
}
