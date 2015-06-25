/* Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include <algorithm>
#include "H5Cpp.h"
#include "hdf5Genome.h"
#include "hdf5DNA.h"
#include "hdf5TopSegment.h"
#include "hdf5BottomSegment.h"
#include "hdf5Sequence.h"
#include "hdf5SequenceIterator.h"
#include "defaultTopSegmentIterator.h"
#include "defaultBottomSegmentIterator.h"
#include "hdf5DNAIterator.h"
#include "defaultColumnIterator.h"
#include "defaultRearrangement.h"
#include "defaultGappedTopSegmentIterator.h"
#include "defaultGappedBottomSegmentIterator.h"

using namespace hal;
using namespace std;
using namespace H5;

const string HDF5Genome::dnaArrayName = "DNA_ARRAY";
const string HDF5Genome::topArrayName = "TOP_ARRAY";
const string HDF5Genome::bottomArrayName = "BOTTOM_ARRAY";
const string HDF5Genome::sequenceIdxArrayName = "SEQIDX_ARRAY";
const string HDF5Genome::sequenceNameArrayName = "SEQNAME_ARRAY";
const string HDF5Genome::metaGroupName = "Meta";
const string HDF5Genome::rupGroupName = "Rup";
const double HDF5Genome::dnaChunkScale = 10.;

HDF5Genome::HDF5Genome(const string& name,
                       HDF5Alignment* alignment,
                       CommonFG* h5Parent,
                       const DSetCreatPropList& dcProps,
                       bool inMemory) :
  _alignment(alignment),
  _h5Parent(h5Parent),
  _name(name),
  _numChildrenInBottomArray(0),
  _totalSequenceLength(0),
  _numChunksInArrayBuffer(inMemory ? 0 : 1),
  _parentCache(NULL)
{
  _dcprops.copy(dcProps);
  assert(!name.empty());
  assert(alignment != NULL && h5Parent != NULL);

  H5::Exception::dontPrint();
  try
  {
    _group = h5Parent->openGroup(name);
    read();
  }
  catch (Exception& e)
  {
    _group = h5Parent->createGroup(name);
  }
  _metaData = new HDF5MetaData(&_group, metaGroupName);
  _rup = new HDF5MetaData(&_group, rupGroupName);

  _totalSequenceLength = _dnaArray.getSize() * 2;
  if (_totalSequenceLength > 0 && _rup->get(rupGroupName) == "1")
  {
    _totalSequenceLength -= 1;
  }
  else if (_totalSequenceLength == 0 && _sequenceIdxArray.getSize() > 0)
  {
    HDF5Sequence lastSeq(this, &_sequenceIdxArray, &_sequenceNameArray,
                          _sequenceNameArray.getSize() - 1);
    _totalSequenceLength = lastSeq.getEndPosition() + 1;
  }

  hsize_t chunk;
  _dcprops.getChunk(1, &chunk);
}


HDF5Genome::~HDF5Genome()
{
  delete _metaData;
  delete _rup;
  deleteSequenceCache();
}

//GENOME INTERFACE

void HDF5Genome::setDimensions(
  const vector<Sequence::Info>& sequenceDimensions,
  bool storeDNAArrays)
{
  _totalSequenceLength = 0;
  hal_size_t totalSeq = sequenceDimensions.size();
  hal_size_t maxName = 0;
  
  // Copy segment dimensions to use the external interface
  vector<Sequence::UpdateInfo> topDimensions;
  topDimensions.reserve(sequenceDimensions.size());
  vector<Sequence::UpdateInfo> bottomDimensions;
  bottomDimensions.reserve(sequenceDimensions.size());

  // Compute summary info from the list of sequence Dimensions
  for (vector<Sequence::Info>::const_iterator i = sequenceDimensions.begin();
       i != sequenceDimensions.end(); 
       ++i)
  {
    _totalSequenceLength += i->_length;
    maxName = max(static_cast<hal_size_t>(i->_name.length()), maxName);
    topDimensions.push_back(
      Sequence::UpdateInfo(i->_name, i->_numTopSegments));
    bottomDimensions.push_back(
      Sequence::UpdateInfo(i->_name, i->_numBottomSegments));
  }

  // Unlink the DNA and segment arrays if they exist (using 
  // exceptions is the only way I know how right now).  Note that
  // the file needs to be refactored to take advantage of the new
  // space. 
  H5::Exception::dontPrint();
  try
  {
    DataSet d = _group.openDataSet(dnaArrayName);
    _group.unlink(dnaArrayName);
  }
  catch (H5::Exception){}
  try
  {
    DataSet d = _group.openDataSet(sequenceIdxArrayName);
    _group.unlink(sequenceIdxArrayName);
  }
  catch (H5::Exception){}
  try
  {
    DataSet d = _group.openDataSet(sequenceNameArrayName);
    _group.unlink(sequenceNameArrayName);
  }
  catch (H5::Exception){}

  if (_totalSequenceLength > 0 && storeDNAArrays == true)
  {
    hal_size_t arrayLength = _totalSequenceLength / 2;
    if (_totalSequenceLength % 2)
    {
      ++arrayLength;
      _rup->set(rupGroupName, "1");
    }
    else
    {
      _rup->set(rupGroupName, "0");
    }
    hsize_t chunk;
    _dcprops.getChunk(1, &chunk);
    // enalarge chunk size because dna bases are so much smaller
    // than segments.  (about 30x). we default to 10x enlargement
    // since the seem to compress about 3x worse.  
    chunk *= dnaChunkScale;
    DSetCreatPropList dnaDC;
    dnaDC.copy(_dcprops);
    dnaDC.setChunk(1, &chunk);
    _dnaArray.create(&_group, dnaArrayName, HDF5DNA::dataType(), 
                     arrayLength, &dnaDC, _numChunksInArrayBuffer);
  }
  if (totalSeq > 0)
  {
    _sequenceIdxArray.create(&_group, sequenceIdxArrayName, 
                             HDF5Sequence::idxDataType(), 
                             totalSeq + 1, &_dcprops, _numChunksInArrayBuffer);

    _sequenceNameArray.create(&_group, sequenceNameArrayName, 
                              HDF5Sequence::nameDataType(maxName + 1), 
                              totalSeq, &_dcprops, _numChunksInArrayBuffer);

    writeSequences(sequenceDimensions);    
  }
  
  // Do the same as above for the segments. 
  setGenomeTopDimensions(topDimensions);
  setGenomeBottomDimensions(bottomDimensions);

  _parentCache = NULL;
  _childCache.clear();
}

void HDF5Genome::updateTopDimensions(
  const vector<Sequence::UpdateInfo>& topDimensions)
{
  loadSequencePosCache();
  loadSequenceNameCache();
  vector<Sequence::UpdateInfo>::const_iterator i;
  map<string, HDF5Sequence*>::iterator cacheIt;
  map<string, const Sequence::UpdateInfo*> inputMap;
  map<string, hal_size_t> currentTopD;
  // copy input into map, checking everything is already present
  for (i = topDimensions.begin(); i != topDimensions.end(); ++i)
  {
    const string& name = i->_name;
    cacheIt = _sequenceNameCache.find(name);
    if (cacheIt == _sequenceNameCache.end())
    {
      throw hal_exception(string("Cannot update sequence ") +
                          name + " because it is not present in "
                          " genome " + getName());
    }
    inputMap.insert(pair<string, const Sequence::UpdateInfo*>(name, &*i));
  }
  // keep a record of the number of segments in each existing 
  // segment (these can get muddled as we add the new ones in the next
  // loop to be sure by getting them in one shot)
  // Note to self: iterating the map in this way skips zero-length 
  // sequences (which are in the separate vector).  This is fine
  // here since we will never update them, but seems like it could be 
  // dangerous if something were to change
  map<hal_size_t, HDF5Sequence*>::iterator posCacheIt;
  map<string, const Sequence::UpdateInfo*>::iterator inputIt;
  for (posCacheIt = _sequencePosCache.begin(); 
       posCacheIt != _sequencePosCache.end(); ++posCacheIt)
  {
    HDF5Sequence* sequence = posCacheIt->second;
    inputIt = inputMap.find(sequence->getName());
    if (inputIt == inputMap.end())
    {
      currentTopD.insert(pair<string, hal_size_t>(
                           sequence->getName(), 
                           sequence->getNumTopSegments()));
    }
  }
  // scan through existing sequences, updating as necessary
  // build summary of all new and unchanged dimensions in newDimensions
  // Note to self: iterating the map in this way skips zero-length 
  // sequences (which are in the separate vector).  This is fine
  // here since we will never update them, but seems like it could be 
  // dangerous if something were to change
  map<string, hal_size_t>::iterator currentIt;
  vector<Sequence::UpdateInfo> newDimensions;
  Sequence::UpdateInfo newInfo;
  hal_size_t topArrayIndex = 0;
  for (posCacheIt = _sequencePosCache.begin(); 
       posCacheIt != _sequencePosCache.end(); ++posCacheIt)
  {
    HDF5Sequence* sequence = posCacheIt->second;
    sequence->setTopSegmentArrayIndex(topArrayIndex);
    inputIt = inputMap.find(sequence->getName());
    if (inputIt != inputMap.end())
    {
      const Sequence::UpdateInfo* updateInfo = inputIt->second;
      newDimensions.push_back(*updateInfo);
    }
    else
    {
      currentIt = currentTopD.find(sequence->getName());
      assert(currentIt != currentTopD.end());
      newInfo._name = posCacheIt->first;
      newInfo._numSegments = currentIt->second;
      newDimensions.push_back(newInfo);
    }
    sequence->setNumTopSegments(newDimensions.back()._numSegments);
    topArrayIndex += newDimensions.back()._numSegments;
  }
  setGenomeTopDimensions(newDimensions);
}

void HDF5Genome::updateBottomDimensions(
  const vector<Sequence::UpdateInfo>& bottomDimensions)
{
  loadSequencePosCache();
  loadSequenceNameCache();
  vector<Sequence::UpdateInfo>::const_iterator i;
  map<string, HDF5Sequence*>::iterator cacheIt;
  map<string, const Sequence::UpdateInfo*> inputMap;
  map<string, hal_size_t> currentBottomD;
  // copy input into map, checking everything is already present
  for (i = bottomDimensions.begin(); i != bottomDimensions.end(); ++i)
  {
    const string& name = i->_name;
    cacheIt = _sequenceNameCache.find(name);
    if (cacheIt == _sequenceNameCache.end())
    {
      throw hal_exception(string("Cannot update sequence ") +
                          name + " because it is not present in "
                          " genome " + getName());
    }
    inputMap.insert(pair<string, const Sequence::UpdateInfo*>(name, &*i));
  }
  // keep a record of the number of segments in each existing 
  // segment (these can get muddled as we add the new ones in the next
  // loop to be sure by getting them in one shot)
  // Note to self: iterating the map in this way skips zero-length 
  // sequences (which are in the separate vector).  This is fine
  // here since we will never update them, but seems like it could be 
  // dangerous if something were to change
  map<hal_size_t, HDF5Sequence*>::iterator posCacheIt;
  map<string, const Sequence::UpdateInfo*>::iterator inputIt;
  for (posCacheIt = _sequencePosCache.begin(); 
       posCacheIt != _sequencePosCache.end(); ++posCacheIt)
  {
    HDF5Sequence* sequence = posCacheIt->second;
    inputIt = inputMap.find(sequence->getName());
    if (inputIt == inputMap.end())
    {
      currentBottomD.insert(pair<string, hal_size_t>(
                           sequence->getName(), 
                           sequence->getNumBottomSegments()));
    }
  }
  // scan through existing sequences, updating as necessary
  // build summary of all new and unchanged dimensions in newDimensions
  // Note to self: iterating the map in this way skips zero-length 
  // sequences (which are in the separate vector).  This is fine
  // here since we will never update them, but seems like it could be 
  // dangerous if something were to change
  map<string, hal_size_t>::iterator currentIt;
  vector<Sequence::UpdateInfo> newDimensions;
  Sequence::UpdateInfo newInfo;
  hal_size_t bottomArrayIndex = 0;
  for (posCacheIt = _sequencePosCache.begin(); 
       posCacheIt != _sequencePosCache.end(); ++posCacheIt)
  {
    HDF5Sequence* sequence = posCacheIt->second;
    sequence->setBottomSegmentArrayIndex(bottomArrayIndex);
    inputIt = inputMap.find(sequence->getName());
    if (inputIt != inputMap.end())
    {
      const Sequence::UpdateInfo* updateInfo = inputIt->second;
      newDimensions.push_back(*updateInfo);
    }
    else
    {
      currentIt = currentBottomD.find(sequence->getName());
      assert(currentIt != currentBottomD.end());
      newInfo._name = posCacheIt->first;
      newInfo._numSegments = currentIt->second;
      newDimensions.push_back(newInfo);
    }
    sequence->setNumBottomSegments(newDimensions.back()._numSegments);
    bottomArrayIndex += newDimensions.back()._numSegments;
  }
  setGenomeBottomDimensions(newDimensions);
}

void HDF5Genome::setGenomeTopDimensions(
  const vector<Sequence::UpdateInfo>& topDimensions)
{
  hal_size_t numTopSegments = 0;
  for (vector<Sequence::UpdateInfo>::const_iterator i = topDimensions.begin();
       i != topDimensions.end(); 
       ++i)
  {
    numTopSegments += i->_numSegments;
  }
  H5::Exception::dontPrint();
  try
  {
    DataSet d = _group.openDataSet(topArrayName);
    _group.unlink(topArrayName);
  }
  catch (H5::Exception){}
  _topArray.create(&_group, topArrayName, HDF5TopSegment::dataType(), 
                   numTopSegments + 1, &_dcprops, _numChunksInArrayBuffer);
  _parentCache = NULL;
}

void HDF5Genome::setGenomeBottomDimensions(
  const vector<Sequence::UpdateInfo>& bottomDimensions)
{
  hal_size_t numBottomSegments = 0;
  for (vector<Sequence::UpdateInfo>::const_iterator i
          = bottomDimensions.begin(); i != bottomDimensions.end(); 
       ++i)
  {
    numBottomSegments += i->_numSegments;
  }
  H5::Exception::dontPrint();
  try
  {
    DataSet d = _group.openDataSet(bottomArrayName);
    _group.unlink(bottomArrayName);
  }
  catch (H5::Exception){}
  hal_size_t numChildren = _alignment->getChildNames(_name).size();
 
  // scale down the chunk size in order to keep chunks proportional to
  // the size of a bottom segment with two children.
  hsize_t chunk;
  _dcprops.getChunk(1, &chunk);  
  double scale = numChildren < 10 ? 1. : 10. / numChildren;
  chunk *= scale;
  DSetCreatPropList botDC;
  botDC.copy(_dcprops);
  botDC.setChunk(1, &chunk);

  _bottomArray.create(&_group, bottomArrayName, 
                      HDF5BottomSegment::dataType(numChildren), 
                      numBottomSegments + 1, &botDC, _numChunksInArrayBuffer);
  _numChildrenInBottomArray = numChildren;
  _childCache.clear();
}

hal_size_t HDF5Genome::getNumSequences() const
{
  assert(_sequenceIdxArray.getSize() == _sequenceNameArray.getSize() + 1);
  return _sequenceNameArray.getSize();
}
   
Sequence* HDF5Genome::getSequence(const string& name)
{
  loadSequenceNameCache();
  Sequence* sequence = NULL;
  map<string, HDF5Sequence*>::iterator mapIt = _sequenceNameCache.find(name);
  if (mapIt != _sequenceNameCache.end())
  {
    sequence = mapIt->second;
  }
  return sequence;
}

const Sequence* HDF5Genome::getSequence(const string& name) const
{
  loadSequenceNameCache();
  const Sequence* sequence = NULL;
  map<string, HDF5Sequence*>::const_iterator mapIt = 
     _sequenceNameCache.find(name);
  if (mapIt != _sequenceNameCache.end())
  {
    sequence = mapIt->second;
  }
  return sequence;
}

Sequence* HDF5Genome::getSequenceBySite(hal_size_t position)
{
  loadSequencePosCache();
  map<hal_size_t, HDF5Sequence*>::iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= (hal_size_t)i->second->getStartPosition())
    {
      assert(position < i->second->getStartPosition() +
             i->second->getSequenceLength());
      return i->second;
    }
  }
  return NULL;
}

const Sequence* HDF5Genome::getSequenceBySite(hal_size_t position) const
{
  loadSequencePosCache();
  map<hal_size_t, HDF5Sequence*>::const_iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= (hal_size_t)i->second->getStartPosition())
    {
      assert(position < i->second->getStartPosition() +
             i->second->getSequenceLength());
      return i->second;
    }
  }
  return NULL;
}

SequenceIteratorPtr HDF5Genome::getSequenceIterator(
  hal_index_t position)
{
  assert(position <= (hal_index_t)_sequenceNameArray.getSize());
  HDF5SequenceIterator* newIt = new HDF5SequenceIterator(this, position);
  return SequenceIteratorPtr(newIt);
}

SequenceIteratorConstPtr HDF5Genome::getSequenceIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)_sequenceNameArray.getSize());
  // genome effectively gets re-consted when returned in the
  // const iterator.  just save doubling up code.
  HDF5SequenceIterator* newIt = new HDF5SequenceIterator(
    const_cast<HDF5Genome*>(this), position);
  return SequenceIteratorConstPtr(newIt);
}

SequenceIteratorConstPtr HDF5Genome::getSequenceEndIterator() const
{
  return getSequenceIterator(getNumSequences());
}

MetaData* HDF5Genome::getMetaData()
{
  return _metaData;
}

const MetaData* HDF5Genome::getMetaData() const
{
  return _metaData;
}

Genome* HDF5Genome::getParent()
{
  if (_parentCache == NULL)
  {
    string parName = _alignment->getParentName(_name);
    if (parName.empty() == false)
    {
      _parentCache = _alignment->openGenome(parName);
    }
  }
  return _parentCache;
}

const Genome* HDF5Genome::getParent() const
{
  if (_parentCache == NULL)
  {
    string parName = _alignment->getParentName(_name);
    if (parName.empty() == false)
    {
      _parentCache = _alignment->openGenome(parName);
    }
  }
  return _parentCache;
}

Genome* HDF5Genome::getChild(hal_size_t childIdx)
{
  assert(childIdx < _numChildrenInBottomArray);
  if (_childCache.size() <= childIdx)
  {
    _childCache.assign(_numChildrenInBottomArray, NULL);
  }
  if (_childCache[childIdx] == NULL)
  {
    vector<string> childNames = _alignment->getChildNames(_name);
    assert(childNames.size() > childIdx);
    _childCache[childIdx] = _alignment->openGenome(childNames.at(childIdx));
  }
  return _childCache[childIdx];
}

const Genome* HDF5Genome::getChild(hal_size_t childIdx) const
{
  assert(childIdx < _numChildrenInBottomArray);
  if (_childCache.size() <= childIdx)
  {
    _childCache.assign(_numChildrenInBottomArray, NULL);
  }
  if (_childCache[childIdx] == NULL)
  {
    vector<string> childNames = _alignment->getChildNames(_name);
    assert(childNames.size() > childIdx);
    _childCache[childIdx] = _alignment->openGenome(childNames.at(childIdx));
  }
  return _childCache[childIdx];
}

hal_size_t HDF5Genome::getNumChildren() const
{
  return _numChildrenInBottomArray;
}

hal_index_t HDF5Genome::getChildIndex(const Genome* child) const
{
  string childName = child->getName();
  vector<string> childNames = _alignment->getChildNames(_name);
  for (hal_size_t i = 0; i < childNames.size(); ++i)
  {
    if (childNames[i] == childName)
    {
      return i;
    }
  }
  return NULL_INDEX;
}

bool HDF5Genome::containsDNAArray() const
{
  return _dnaArray.getSize() > 0;
}

const Alignment* HDF5Genome::getAlignment() const
{
  return _alignment;
}

// SEGMENTED SEQUENCE INTERFACE

const string& HDF5Genome::getName() const
{
  return _name;
}

hal_size_t HDF5Genome::getSequenceLength() const
{
  return _totalSequenceLength;
}

hal_size_t HDF5Genome::getNumTopSegments() const
{
  hal_size_t arraySize = _topArray.getSize();
  return arraySize > 0 ? arraySize - 1 : 0;
}

hal_size_t HDF5Genome::getNumBottomSegments() const
{
  hal_size_t arraySize = _bottomArray.getSize();
  return arraySize > 0 ? arraySize - 1 : 0;
}

TopSegmentIteratorPtr HDF5Genome::getTopSegmentIterator(hal_index_t position)
{
  assert(position <= (hal_index_t)getNumTopSegments());
  HDF5TopSegment* newSeg = new HDF5TopSegment(this, &_topArray, position);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultTopSegmentIterator* newIt = new DefaultTopSegmentIterator(newSeg);
  return TopSegmentIteratorPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Genome::getTopSegmentIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)getNumTopSegments());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  HDF5TopSegment* newSeg = new HDF5TopSegment(genome, &genome->_topArray, 
                                              position);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultTopSegmentIterator* newIt = new DefaultTopSegmentIterator(newSeg);
  return TopSegmentIteratorConstPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Genome::getTopSegmentEndIterator() const
{
  return getTopSegmentIterator(getNumTopSegments());
}

BottomSegmentIteratorPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position)
{
  assert(position <= (hal_index_t)getNumBottomSegments());
  HDF5BottomSegment* newSeg = new HDF5BottomSegment(this, &_bottomArray, 
                                                    position);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultBottomSegmentIterator* newIt = new DefaultBottomSegmentIterator(newSeg);
  return BottomSegmentIteratorPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)getNumBottomSegments());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  HDF5BottomSegment* newSeg = new HDF5BottomSegment(genome, &genome->_bottomArray, 
                                                    position);
  // ownership of newSeg is passed into newIt, whose lifespan is 
  // governed by the returned smart pointer
  DefaultBottomSegmentIterator* newIt = new DefaultBottomSegmentIterator(newSeg);
  return BottomSegmentIteratorConstPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Genome::getBottomSegmentEndIterator() const
{
  return getBottomSegmentIterator(getNumBottomSegments());
}
   
DNAIteratorPtr HDF5Genome::getDNAIterator(hal_index_t position)
{
  assert(position / 2 <= (hal_index_t)_dnaArray.getSize());
  HDF5DNAIterator* newIt = new HDF5DNAIterator(this, position);
  return DNAIteratorPtr(newIt);
}

DNAIteratorConstPtr HDF5Genome::getDNAIterator(hal_index_t position) const
{
  assert(_dnaArray.getSize() == 0 || 
         position / 2 <= (hal_index_t)_dnaArray.getSize());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  const HDF5DNAIterator* newIt = new HDF5DNAIterator(genome, position);
  return DNAIteratorConstPtr(newIt);
}

DNAIteratorConstPtr HDF5Genome::getDNAEndIterator() const
{
  return getDNAIterator(getSequenceLength());
}

ColumnIteratorConstPtr HDF5Genome::getColumnIterator(
  const set<const Genome*>* targets, hal_size_t maxInsertLength, 
  hal_index_t position, hal_index_t lastPosition, bool noDupes,
  bool noAncestors, bool reverseStrand, bool unique, bool onlyOrthologs) const
{
  hal_index_t lastIdx = lastPosition;
  if (lastPosition == NULL_INDEX)
  {
    lastIdx = (hal_index_t)(getSequenceLength() - 1);
  }
  if (position < 0 || 
      lastPosition >= (hal_index_t)(getSequenceLength()))
  {
    stringstream ss;
    ss << "HDF5Genome::getColumnIterator: input indices "
       << "(" << position << ", " << lastPosition << ") out of bounds";
    throw hal_exception(ss.str());
  }
  const DefaultColumnIterator* newIt = 
     new DefaultColumnIterator(this, targets, position, lastIdx, 
                               maxInsertLength, noDupes, noAncestors,
                               reverseStrand, unique, onlyOrthologs);
  return ColumnIteratorConstPtr(newIt);
}

void HDF5Genome::getString(string& outString) const
{
  getSubString(outString, 0, getSequenceLength());
}

void HDF5Genome::setString(const string& inString)
{
  setSubString(inString, 0, getSequenceLength());
}

void HDF5Genome::getSubString(string& outString, hal_size_t start,
                              hal_size_t length) const
{
  outString.resize(length);
  HDF5DNAIterator dnaIt(const_cast<HDF5Genome*>(this), start);
  dnaIt.readString(outString, length);
}

void HDF5Genome::setSubString(const string& inString, 
                              hal_size_t start,
                              hal_size_t length)
{
  if (length != inString.length())
  {
    throw hal_exception(string("setString: input string has differnt") +
                               "length from target string in genome");
  }
  HDF5DNAIterator dnaIt(this, start);
  dnaIt.writeString(inString, length);
}

RearrangementPtr HDF5Genome::getRearrangement(hal_index_t position,
                                              hal_size_t gapLengthThreshold,
                                              double nThreshold,
                                              bool atomic) const
{
  assert(position >= 0 && position < (hal_index_t)getNumTopSegments());
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(position);  
  DefaultRearrangement* rea = new DefaultRearrangement(this,
                                                       gapLengthThreshold,
                                                       nThreshold,
                                                       atomic);
  rea->identifyFromLeftBreakpoint(top);
  return RearrangementPtr(rea);
}

GappedTopSegmentIteratorConstPtr HDF5Genome::getGappedTopSegmentIterator(
  hal_index_t i, hal_size_t gapThreshold, bool atomic) const
{
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(i);  
  DefaultGappedTopSegmentIterator* gt = 
     new DefaultGappedTopSegmentIterator(top, gapThreshold, atomic);
  return GappedTopSegmentIteratorConstPtr(gt);
}

GappedBottomSegmentIteratorConstPtr HDF5Genome::getGappedBottomSegmentIterator(
  hal_index_t i, hal_size_t childIdx, hal_size_t gapThreshold,
  bool atomic) const
{
  BottomSegmentIteratorConstPtr bot = getBottomSegmentIterator(i);  
  DefaultGappedBottomSegmentIterator* gb = 
     new DefaultGappedBottomSegmentIterator(bot, childIdx, gapThreshold, 
                                            atomic);
  return GappedBottomSegmentIteratorConstPtr(gb);
}

  
// LOCAL NON-INTERFACE METHODS

void HDF5Genome::write()
{
  _dnaArray.write();
  _topArray.write();
  _bottomArray.write();
  _metaData->write();
  _rup->write();
  _sequenceIdxArray.write();
  _sequenceNameArray.write();
}

void HDF5Genome::read()
{
  H5::Exception::dontPrint();
  try
  {
    _group.openDataSet(dnaArrayName);
    _dnaArray.load(&_group, dnaArrayName, _numChunksInArrayBuffer);    
  }
  catch (H5::Exception){}

  try
  {
    _group.openDataSet(topArrayName);
    _topArray.load(&_group, topArrayName, _numChunksInArrayBuffer);
  }
  catch (H5::Exception){}
  try
  {
    _group.openDataSet(bottomArrayName);
    _bottomArray.load(&_group, bottomArrayName, _numChunksInArrayBuffer);
    _numChildrenInBottomArray = 
       HDF5BottomSegment::numChildrenFromDataType(_bottomArray.getDataType());
  }
  catch (H5::Exception){}

  deleteSequenceCache();
  try
  {
    _group.openDataSet(sequenceIdxArrayName);
    _sequenceIdxArray.load(&_group, sequenceIdxArrayName, 
                           _numChunksInArrayBuffer);
  }
  catch (H5::Exception){}
  try
  {
    _group.openDataSet(sequenceNameArrayName);
    _sequenceNameArray.load(&_group, sequenceNameArrayName, 
                            _numChunksInArrayBuffer);
  }
  catch (H5::Exception){}

  readSequences();
}

void HDF5Genome::readSequences()
{
  deleteSequenceCache();
}

void HDF5Genome::deleteSequenceCache()
{
  if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0)
  {
    map<hal_size_t, HDF5Sequence*>::iterator i;
    for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i)
    {
      delete i->second;
    }
    vector<HDF5Sequence*>::iterator z;
    for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z)
    {
      delete *z;
    }
  }
  else if (_sequenceNameCache.size() > 0)
  {
    map<string, HDF5Sequence*>::iterator i;
    for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i)
    {
      delete i->second;
    }
  }
  _sequencePosCache.clear();
  _zeroLenPosCache.clear();
  _sequenceNameCache.clear(); // I share my pointers with above. 
}

void HDF5Genome::loadSequencePosCache() const
{
  if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0)
  {
    return;
  }
  hal_size_t totalReadLen = 0;
  hal_size_t numSequences = _sequenceNameArray.getSize();
  
  if (_sequenceNameCache.size() > 0)
  {
    assert(_sequenceNameCache.size() == numSequences);
    map<std::string, HDF5Sequence*>::iterator i;
    for (i = _sequenceNameCache.begin(); i != _sequenceNameCache.end(); ++i)
    {
      if (i->second->getSequenceLength() > 0)
      {
        _sequencePosCache.insert(pair<hal_size_t, HDF5Sequence*>(
                                   i->second->getStartPosition() +
                                   i->second->getSequenceLength(), i->second));
        totalReadLen += i->second->getSequenceLength();
      }
      else
      {
        _zeroLenPosCache.push_back(i->second);
      }
    }
  }
  else
  {
    for (hal_size_t i = 0; i < numSequences; ++i)
    {
      HDF5Sequence* seq = 
         new HDF5Sequence(const_cast<HDF5Genome*>(this),
                          const_cast<HDF5ExternalArray*>(&_sequenceIdxArray),
                          const_cast<HDF5ExternalArray*>(&_sequenceNameArray),
                          i);
      if (seq->getSequenceLength() > 0)
      {
        _sequencePosCache.insert(
          pair<hal_size_t, HDF5Sequence*>(seq->getStartPosition() +
                                          seq->getSequenceLength(), seq));
        totalReadLen += seq->getSequenceLength();
      }
      else
      {
        _zeroLenPosCache.push_back(seq);
      }
    }
  }
  if (_totalSequenceLength > 0 && totalReadLen != _totalSequenceLength)
  {
    stringstream ss;
    ss << "Sequences for genome " << getName() << " have total length " 
       << totalReadLen << " but the (non-zero) DNA array contains "
       << _totalSequenceLength << " elements. This is an internal error "
       << "or the file is corrupt.";
    throw hal_exception(ss.str());
  }
}

void HDF5Genome::loadSequenceNameCache() const
{
  if (_sequenceNameCache.size() > 0)
  {
    return;
  }
  hal_size_t numSequences = _sequenceNameArray.getSize();
  
  if (_sequencePosCache.size() > 0 || _zeroLenPosCache.size() > 0)
  {
    assert(_sequencePosCache.size() + _zeroLenPosCache.size() == numSequences);
    map<hal_size_t, HDF5Sequence*>::iterator i;
    for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i)
    {
      _sequenceNameCache.insert(pair<string, HDF5Sequence*>(
                                  i->second->getName(), i->second));
    }
    vector<HDF5Sequence*>::iterator z;
    for (z = _zeroLenPosCache.begin(); z != _zeroLenPosCache.end(); ++z)
    {
      _sequenceNameCache.insert(pair<string, HDF5Sequence*>(
                                  (*z)->getName(), (*z)));
    }
  }
  else
  {
    for (hal_size_t i = 0; i < numSequences; ++i)
    {
      HDF5Sequence* seq = 
         new HDF5Sequence(const_cast<HDF5Genome*>(this),
                          const_cast<HDF5ExternalArray*>(&_sequenceIdxArray),
                          const_cast<HDF5ExternalArray*>(&_sequenceNameArray),
                          i);

      _sequenceNameCache.insert(
        pair<string, HDF5Sequence*>(seq->getName(), seq));
    }
  }
}
  
void HDF5Genome::writeSequences(const vector<Sequence::Info>&
                                sequenceDimensions)
{
  deleteSequenceCache();
  vector<Sequence::Info>::const_iterator i;
  hal_size_t startPosition = 0;
  hal_size_t topArrayIndex = 0;
  hal_size_t bottomArrayIndex = 0;
  for (i = sequenceDimensions.begin(); i != sequenceDimensions.end(); ++i)
  {
    // Copy segment into HDF5 array
    HDF5Sequence* seq = new HDF5Sequence(this, &_sequenceIdxArray,
                                         &_sequenceNameArray,
                                         i - sequenceDimensions.begin());
    // write all the Sequence::Info into the hdf5 sequence record
    seq->set(startPosition, *i, topArrayIndex, bottomArrayIndex);
    // Keep the object pointer in our caches
    if (seq->getSequenceLength() > 0)
    {
      _sequencePosCache.insert(
        pair<hal_size_t, HDF5Sequence*>(startPosition + i->_length, seq));
    }
    else
    {
      _zeroLenPosCache.push_back(seq);
    }
    _sequenceNameCache.insert(pair<string, HDF5Sequence*>(i->_name, seq));
    startPosition += i->_length;
    topArrayIndex += i->_numTopSegments;
    bottomArrayIndex += i->_numBottomSegments;
  }  
}

void HDF5Genome::resetBranchCaches()
{
  _parentCache = NULL;
  _childCache.clear();
}
