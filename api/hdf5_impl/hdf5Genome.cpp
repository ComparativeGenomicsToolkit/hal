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
#include "hdf5TopSegmentIterator.h"
#include "hdf5BottomSegmentIterator.h"
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
const string HDF5Genome::sequenceArrayName = "SEQUENCE_ARRAY";
const string HDF5Genome::metaGroupName = "Meta";

HDF5Genome::HDF5Genome(const string& name,
                       HDF5Alignment* alignment,
                       CommonFG* h5Parent,
                       const DSetCreatPropList& dcProps) :
  _alignment(alignment),
  _h5Parent(h5Parent),
  _name(name),
  _dcprops(dcProps),
  _numChildrenInBottomArray(0),
  _parentCache(NULL)
{
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
}


HDF5Genome::~HDF5Genome()
{
  delete _metaData;
  deleteSequenceCache();
}

//GENOME INTERFACE

void HDF5Genome::setDimensions(
  const vector<Sequence::Info>& sequenceDimensions)
{
  hal_size_t totalSequenceLength = 0;
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
    totalSequenceLength += i->_length;
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
    DataSet d = _group.openDataSet(sequenceArrayName);
    _group.unlink(sequenceArrayName);
  }
  catch (H5::Exception){}
  if (totalSequenceLength > 0)
  {
    _dnaArray.create(&_group, dnaArrayName, HDF5DNA::dataType(), 
                     totalSequenceLength, _dcprops);
  }
  if (totalSeq > 0)
  {
    _sequenceArray.create(&_group, sequenceArrayName, 
                          // pad names a bit to allow renaming
                          HDF5Sequence::dataType(maxName + 32), 
                          totalSeq, DSetCreatPropList());
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
  vector<Sequence::UpdateInfo>::const_iterator i;
  map<string, HDF5Sequence*>::iterator cacheIt;
  map<string, const Sequence::UpdateInfo*> inputMap;
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
  // scan through existing sequences, updating as necessary
  // build summary of all new and unchanged dimensions in newDimensions
  map<hal_size_t, HDF5Sequence*>::iterator posCacheIt;
  map<string, const Sequence::UpdateInfo*>::iterator inputIt;
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
      sequence->setNumTopSegments(updateInfo->_numSegments);
      newDimensions.push_back(*updateInfo);
    }
    else
    {
      newInfo._name = posCacheIt->first;
      newInfo._numSegments = sequence->getNumTopSegments();
      newDimensions.push_back(newInfo);
    }
    topArrayIndex += newDimensions.back()._numSegments;
  }
  setGenomeTopDimensions(newDimensions);
}

void HDF5Genome::updateBottomDimensions(
  const vector<Sequence::UpdateInfo>& bottomDimensions)
{
  vector<Sequence::UpdateInfo>::const_iterator i;
  map<string, HDF5Sequence*>::iterator cacheIt;
  map<string, const Sequence::UpdateInfo*> inputMap;
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
  // scan through existing sequences, updating as necessary
  // build summary of all new and unchanged dimensions in newDimensions
  map<hal_size_t, HDF5Sequence*>::iterator posCacheIt;
  map<string, const Sequence::UpdateInfo*>::iterator inputIt;
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
      sequence->setNumBottomSegments(updateInfo->_numSegments);
      newDimensions.push_back(*updateInfo);
    }
    else
    {
      newInfo._name = posCacheIt->first;
      newInfo._numSegments = sequence->getNumBottomSegments();
      newDimensions.push_back(newInfo);
    }
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
  if (numTopSegments > 0)
  {
    _topArray.create(&_group, topArrayName, HDF5TopSegment::dataType(), 
                     numTopSegments, _dcprops);
  }
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
  if (numBottomSegments > 0)
  {
    _bottomArray.create(&_group, bottomArrayName, 
                        HDF5BottomSegment::dataType(numChildren), 
                        numBottomSegments, _dcprops);
    _numChildrenInBottomArray = numChildren;
  }
  _childCache.clear();
}

hal_size_t HDF5Genome::getNumSequences() const
{
  return _sequenceArray.getSize();
}
   
Sequence* HDF5Genome::getSequence(const string& name)
{
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
  map<hal_size_t, HDF5Sequence*>::iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= i->second->getStartPosition())
    {
      return i->second;
    }
  }
  return NULL;
}

const Sequence* HDF5Genome::getSequenceBySite(hal_size_t position) const
{
  map<hal_size_t, HDF5Sequence*>::const_iterator i;
  i = _sequencePosCache.upper_bound(position);
  if (i != _sequencePosCache.end())
  {
    if (position >= i->second->getStartPosition())
    {
      return i->second;
    }
  }
  return NULL;
}

SequenceIteratorPtr HDF5Genome::getSequenceIterator(
  hal_index_t position)
{
  assert(position <= (hal_index_t)_sequenceArray.getSize());
  HDF5SequenceIterator* newIt = new HDF5SequenceIterator(this, position);
  return SequenceIteratorPtr(newIt);
}

SequenceIteratorConstPtr HDF5Genome::getSequenceIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)_sequenceArray.getSize());
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
  if (_childCache.size() <= childIdx || 
      _childCache[childIdx] == NULL)
  {
    vector<string> childNames = _alignment->getChildNames(_name);
    _childCache.resize(childNames.size());
    for (size_t i = 0; i < _childCache.size(); ++i)
    {
      _childCache[i] = _alignment->openGenome(childNames.at(i));
    }
  }
  return _childCache[childIdx];
}

const Genome* HDF5Genome::getChild(hal_size_t childIdx) const
{
  assert(childIdx < _numChildrenInBottomArray);
  if (_childCache.size() <= childIdx || 
      _childCache[childIdx] == NULL)
  {
    vector<string> childNames = _alignment->getChildNames(_name);
    _childCache.resize(childNames.size());
    for (size_t i = 0; i < _childCache.size(); ++i)
    {
      _childCache[i] = _alignment->openGenome(childNames.at(i));
    }
  }
  return _childCache[childIdx];
}

hal_size_t HDF5Genome::getNumChildren() const
{
  return _numChildrenInBottomArray;
}

hal_index_t HDF5Genome::getChildIndex(const Genome* child) const
{
  for (hal_size_t i = 0; i < _numChildrenInBottomArray; ++i)
  {
    if (getChild(i) == child)
    {
      return i;
    }
  }
  return NULL_INDEX;
}

// SEGMENTED SEQUENCE INTERFACE

const string& HDF5Genome::getName() const
{
  return _name;
}

hal_size_t HDF5Genome::getSequenceLength() const
{
  return _dnaArray.getSize();
}

hal_size_t HDF5Genome::getNumTopSegments() const
{
  return _topArray.getSize();
}

hal_size_t HDF5Genome::getNumBottomSegments() const
{
  return _bottomArray.getSize();
}

TopSegmentIteratorPtr HDF5Genome::getTopSegmentIterator(hal_index_t position)
{
  assert(position <= (hal_index_t)_topArray.getSize());
  HDF5TopSegmentIterator* newIt = new HDF5TopSegmentIterator(this, position);                                                             
  return TopSegmentIteratorPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Genome::getTopSegmentIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)_topArray.getSize());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  const HDF5TopSegmentIterator* newIt = 
     new HDF5TopSegmentIterator(genome, position);
  return TopSegmentIteratorConstPtr(newIt);
}

TopSegmentIteratorConstPtr HDF5Genome::getTopSegmentEndIterator() const
{
  return getTopSegmentIterator(getNumTopSegments());
}

BottomSegmentIteratorPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position)
{
  assert(position <= (hal_index_t)_bottomArray.getSize());
  HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(this, position);
  return BottomSegmentIteratorPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position) const
{
  assert(position <= (hal_index_t)_bottomArray.getSize());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  const HDF5BottomSegmentIterator* newIt = 
     new HDF5BottomSegmentIterator(genome, position);
  return BottomSegmentIteratorConstPtr(newIt);
}

BottomSegmentIteratorConstPtr HDF5Genome::getBottomSegmentEndIterator() const
{
  return getBottomSegmentIterator(getNumBottomSegments());
}
   
DNAIteratorPtr HDF5Genome::getDNAIterator(hal_index_t position)
{
  assert(position <= (hal_index_t)_dnaArray.getSize());
  HDF5DNAIterator* newIt = new HDF5DNAIterator(this, position);
  return DNAIteratorPtr(newIt);
}

DNAIteratorConstPtr HDF5Genome::getDNAIterator(hal_index_t position) const
{
  assert(position <= (hal_index_t)_dnaArray.getSize());
  HDF5Genome* genome = const_cast<HDF5Genome*>(this);
  const HDF5DNAIterator* newIt = new HDF5DNAIterator(genome, position);
  return DNAIteratorConstPtr(newIt);
}

DNAIteratorConstPtr HDF5Genome::getDNAEndIterator() const
{
  return getDNAIterator(getSequenceLength());
}

ColumnIteratorConstPtr HDF5Genome::getColumnIterator(
  const Genome* root, hal_size_t maxInsertLength, hal_index_t position,
  hal_index_t lastPosition, bool noDupes) const
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
     new DefaultColumnIterator(this, root, position, lastIdx, 
                               maxInsertLength, noDupes);
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

RearrangementPtr HDF5Genome::getRearrangement(hal_index_t position) const
{
  assert(position >= 0 && position < (hal_index_t)getNumTopSegments());
  TopSegmentIteratorConstPtr top = getTopSegmentIterator(position);  
  DefaultRearrangement* rea = new DefaultRearrangement(this);
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
  _sequenceArray.write();
}

void HDF5Genome::read()
{
  H5::Exception::dontPrint();
  try
  {
    _group.openDataSet(dnaArrayName);
    _dnaArray.load(&_group, dnaArrayName);
  }
  catch (H5::Exception){}
  try
  {
    _group.openDataSet(topArrayName);
    _topArray.load(&_group, topArrayName);
  }
  catch (H5::Exception){}
  try
  {
    _group.openDataSet(bottomArrayName);
    _bottomArray.load(&_group, bottomArrayName);
    _numChildrenInBottomArray = 
       HDF5BottomSegment::numChildrenFromDataType(_bottomArray.getDataType());
  }
  catch (H5::Exception){}
  deleteSequenceCache();
  try
  {
    _group.openDataSet(sequenceArrayName);
    _sequenceArray.load(&_group, sequenceArrayName);
    readSequences();
  }
  catch (H5::Exception){}
}

void HDF5Genome::readSequences()
{
  deleteSequenceCache();
  hal_size_t numSequences = _sequenceArray.getSize();
  for (hal_size_t i = 0; i < numSequences; ++i)
  {
    HDF5Sequence* seq = new HDF5Sequence(this, &_sequenceArray, i);
    _sequencePosCache.insert(
      pair<hal_size_t, HDF5Sequence*>(seq->getStartPosition() +
                                      seq->getSequenceLength(), seq));
    _sequenceNameCache.insert(
      pair<string, HDF5Sequence*>(seq->getName(), seq));

  }
}

void HDF5Genome::deleteSequenceCache()
{
  map<hal_size_t, HDF5Sequence*>::iterator i;
  for (i = _sequencePosCache.begin(); i != _sequencePosCache.end(); ++i)
  {
    delete i->second;
  }
  _sequencePosCache.clear();
  _sequenceNameCache.clear(); // I share my pointers with above. 
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
    HDF5Sequence* seq = new HDF5Sequence(this, &_sequenceArray, 
                                         i - sequenceDimensions.begin());
    // write all the Sequence::Info into the hdf5 sequence record
    seq->set(startPosition, *i, topArrayIndex, bottomArrayIndex);
    // Keep the object pointer in our caches
    _sequencePosCache.insert(
      pair<hal_size_t, HDF5Sequence*>(startPosition + i->_length, seq));
    _sequenceNameCache.insert(pair<string, HDF5Sequence*>(i->_name, seq));
    startPosition += i->_length;
    topArrayIndex += i->_numTopSegments;
    bottomArrayIndex += i->_numBottomSegments;
  }  
}
