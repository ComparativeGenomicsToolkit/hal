/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <cstring>
#include <sstream>
#include <set>
#include <iostream>
#include "hdf5Sequence.h"
#include "hdf5DNAIterator.h"
#include "defaultTopSegmentIterator.h"
#include "defaultBottomSegmentIterator.h"
#include "defaultColumnIterator.h"
#include "defaultRearrangement.h"
#include "defaultGappedTopSegmentIterator.h"
#include "defaultGappedBottomSegmentIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

const size_t HDF5Sequence::startOffset = 0;
const size_t HDF5Sequence::topSegmentArrayIndexOffset = sizeof(hal_size_t);
const size_t HDF5Sequence::bottomSegmentArrayIndexOffset = topSegmentArrayIndexOffset + sizeof(hal_size_t);
const size_t HDF5Sequence::totalSize = bottomSegmentArrayIndexOffset + sizeof(hal_size_t);

HDF5Sequence::HDF5Sequence(HDF5Genome* genome,
                           HDF5ExternalArray* idxArray,
                           HDF5ExternalArray* nameArray,
                           hal_index_t index) :
  _idxArray(idxArray),
  _nameArray(nameArray),
  _index(index),
  _genome(genome)
{

}

HDF5Sequence::~HDF5Sequence()
{
  
}

H5::CompType HDF5Sequence::idxDataType()
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  CompType dataType(totalSize);
  dataType.insertMember("start", startOffset, PredType::NATIVE_HSIZE);
  dataType.insertMember("topSegmentArrayIndexOffset", 
                        topSegmentArrayIndexOffset, 
                        PredType::NATIVE_HSIZE);
  dataType.insertMember("bottomSegmentArrayIndexOffset",
                        bottomSegmentArrayIndexOffset, 
                        PredType::NATIVE_HSIZE);
  return dataType;
}

H5::StrType HDF5Sequence::nameDataType(hal_size_t maxNameLength)
{
  // the in-memory representations and hdf5 representations 
  // don't necessarily have to be the same, but it simplifies 
  // testing for now. 
  assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
  assert(PredType::NATIVE_CHAR.getSize() == sizeof(char));
  StrType strType(PredType::NATIVE_CHAR, (maxNameLength + 1) * sizeof(char));
//  CompType dataType(nameOffset + strType.getSize());
//  dataType.insertMember("name", nameOffset, strType);
  return strType;
//  return dataType;
}

// SEQUENCE INTERFACE
string HDF5Sequence::getName() const
{
  return _nameArray->get(_index);
}

string HDF5Sequence::getFullName() const
{
  assert(_genome != NULL);
  return _genome->getName() + '.' + getName();
}

const Genome* HDF5Sequence::getGenome() const
{
  return _genome;
}

Genome* HDF5Sequence::getGenome()
{
  return _genome;
}

hal_index_t HDF5Sequence::getStartPosition() const
{
  return _idxArray->getValue<hal_size_t>(_index, startOffset);
}

hal_index_t HDF5Sequence::getEndPosition() const
{
  return _idxArray->getValue<hal_size_t>(_index + 1 , startOffset) - 1;
}

hal_index_t HDF5Sequence::getArrayIndex() const
{
  return _index;
}

hal_index_t HDF5Sequence::getTopSegmentArrayIndex() const
{
  return (hal_index_t)
     _idxArray->getValue<hal_size_t>(_index, topSegmentArrayIndexOffset);
}

hal_index_t HDF5Sequence::getBottomSegmentArrayIndex() const
{
  return (hal_index_t)
     _idxArray->getValue<hal_size_t>(_index, bottomSegmentArrayIndexOffset);
}

// SEGMENTED SEQUENCE INTERFACE

hal_size_t HDF5Sequence::getSequenceLength() const
{
  hal_size_t len = _idxArray->getValue<hal_size_t>(_index, startOffset);
  hal_size_t nlen = _idxArray->getValue<hal_size_t>(_index + 1, startOffset);
  assert(nlen >= len);
  return (nlen - len);
}

hal_size_t HDF5Sequence::getNumTopSegments() const
{
  hal_index_t idx = 
     _idxArray->getValue<hal_size_t>(_index, topSegmentArrayIndexOffset);
  hal_index_t nextIdx = 
     _idxArray->getValue<hal_size_t>(_index + 1, topSegmentArrayIndexOffset);
  assert(nextIdx >= idx);
  return nextIdx - idx;
}

hal_size_t HDF5Sequence::getNumBottomSegments() const
{
  hal_index_t idx = 
     _idxArray->getValue<hal_size_t>(_index, bottomSegmentArrayIndexOffset);
  hal_index_t nextIdx = 
     _idxArray->getValue<hal_size_t>(_index + 1, bottomSegmentArrayIndexOffset);
  assert(nextIdx >= idx);
  return nextIdx - idx;
}

TopSegmentIteratorPtr HDF5Sequence::getTopSegmentIterator(
  hal_index_t position)
{
  hal_size_t idx = position + getTopSegmentArrayIndex();
  return _genome->getTopSegmentIterator(idx);
}

TopSegmentIteratorConstPtr HDF5Sequence::getTopSegmentIterator(
  hal_index_t position) const
{
  hal_size_t idx = position + getTopSegmentArrayIndex();
  return _genome->getTopSegmentIterator(idx);
}

TopSegmentIteratorConstPtr HDF5Sequence::getTopSegmentEndIterator() const
{
  return getTopSegmentIterator(getNumTopSegments());
}

BottomSegmentIteratorPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position)
{
  hal_size_t idx = position + getBottomSegmentArrayIndex();
  return _genome->getBottomSegmentIterator(idx);
}

BottomSegmentIteratorConstPtr HDF5Sequence::getBottomSegmentIterator(
  hal_index_t position) const
{
  hal_size_t idx = position + getBottomSegmentArrayIndex();
  return _genome->getBottomSegmentIterator(idx);
}

BottomSegmentIteratorConstPtr HDF5Sequence::getBottomSegmentEndIterator() const
{
  return getBottomSegmentIterator(getNumBottomSegments());
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
  return getDNAIterator(getSequenceLength());
}

ColumnIteratorConstPtr HDF5Sequence::getColumnIterator(
  const std::set<const Genome*>* targets, hal_size_t maxInsertLength, 
  hal_index_t position, hal_index_t lastPosition, bool noDupes,
  bool noAncestors, bool reverseStrand, bool unique, bool onlyOrthologs) const
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
    ss << "HDF5Sequence::getColumnIterators: input indices "
       << "(" << position << ", " << lastPosition << ") out of bounds";
    throw hal_exception(ss.str());
  }
  const DefaultColumnIterator* newIt = 
     new DefaultColumnIterator(getGenome(), targets, idx, lastIdx, 
                               maxInsertLength, noDupes, noAncestors,
                               reverseStrand, unique, onlyOrthologs);
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

RearrangementPtr HDF5Sequence::getRearrangement(hal_index_t position,
                                                hal_size_t gapLengthThreshold,
                                                double nThreshold,
                                                bool atomic) const
{
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(position);  
  DefaultRearrangement* rea = new DefaultRearrangement(getGenome(),
                                                       gapLengthThreshold,
                                                       nThreshold,
                                                       atomic);
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
  _idxArray->setValue(_index, startOffset, startPosition);
  _idxArray->setValue(_index + 1, startOffset, 
                      startPosition + sequenceInfo._length);
  _idxArray->setValue(_index, topSegmentArrayIndexOffset, topSegmentStartIndex);
  _idxArray->setValue(_index, bottomSegmentArrayIndexOffset, 
                      bottomSegmentStartIndex);
  _idxArray->setValue(_index + 1, topSegmentArrayIndexOffset, 
                      topSegmentStartIndex + sequenceInfo._numTopSegments);
  _idxArray->setValue(_index + 1, bottomSegmentArrayIndexOffset,
                      bottomSegmentStartIndex + sequenceInfo._numBottomSegments);
  char* arrayBuffer = _nameArray->getUpdate(_index);
  strcpy(arrayBuffer, sequenceInfo._name.c_str());

  assert(getStartPosition() == (hal_index_t)startPosition);
  assert(getNumTopSegments() == sequenceInfo._numTopSegments);
  assert(getNumBottomSegments() == sequenceInfo._numBottomSegments);
  assert(getSequenceLength() == sequenceInfo._length);
}

// These functions look dangerous.  Dont think they're used.
void HDF5Sequence::setNumTopSegments(hal_size_t numTopSegments)
{
  _idxArray->setValue(_index + 1, topSegmentArrayIndexOffset, 
                      getTopSegmentArrayIndex() + numTopSegments);
}

void HDF5Sequence::setNumBottomSegments(hal_size_t numBottomSegments)
{
  _idxArray->setValue(_index + 1, bottomSegmentArrayIndexOffset, 
                      getBottomSegmentArrayIndex() + numBottomSegments);
}

void HDF5Sequence::setTopSegmentArrayIndex(hal_size_t topIndex)
{
  _idxArray->setValue(_index, topSegmentArrayIndexOffset, topIndex);
}

void HDF5Sequence::setBottomSegmentArrayIndex(hal_size_t bottomIndex)
{
  _idxArray->setValue(_index, bottomSegmentArrayIndexOffset, bottomIndex);
}
