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
#include "hdf5DNAIterator.h"

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
                          totalSeq, _dcprops);
    writeSequences(sequenceDimensions);
    
  }
  
  // Do the same as above for the segments. 
  setTopDimensions(topDimensions);
  setBottomDimensions(bottomDimensions);

  _parentCache = NULL;
  _childCache.clear();
}

void HDF5Genome::setTopDimensions(
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

void HDF5Genome::setBottomDimensions(
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
  map<string, HDF5Sequence*>::iterator mapIt = _sequenceNameCache.find(name);
  if (mapIt == _sequenceNameCache.end())
  {
    throw hal_exception(name + ":  sequence not in genome");
  }
  return mapIt->second;
}

const Sequence* HDF5Genome::getSequence(const string& name) const
{
  map<string, HDF5Sequence*>::const_iterator mapIt = 
     _sequenceNameCache.find(name);
  if (mapIt == _sequenceNameCache.end())
  {
    throw hal_exception(name + ":  sequence not in genome");
  }
  return mapIt->second;
}

Sequence* HDF5Genome::getSequenceBySite(hal_index_t position)
{
  return NULL;
}

const Sequence* HDF5Genome::getSequenceBySite(hal_index_t position) const
{
  return NULL;
}

SequenceIteratorPtr HDF5Genome::getSequenceIterator(
  hal_index_t position)
{

}

SequenceIteratorConstPtr HDF5Genome::getSequenceIterator(
  hal_index_t position) const
{

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
    _parentCache = _alignment->openGenome(
      _alignment->getParentName(_name));
  }
  return _parentCache;
}

const Genome* HDF5Genome::getParent() const
{
  if (_parentCache == NULL)
  {
    _parentCache = _alignment->openGenome(
      _alignment->getParentName(_name));
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
  return TopSegmentIteratorPtr(0);
}

TopSegmentIteratorConstPtr HDF5Genome::getTopSegmentIterator(
  hal_index_t position) const
{
  return TopSegmentIteratorConstPtr(0);
}

BottomSegmentIteratorPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position)
{
  return BottomSegmentIteratorPtr(0);
}

BottomSegmentIteratorConstPtr HDF5Genome::getBottomSegmentIterator(
  hal_index_t position) const
{
  return BottomSegmentIteratorConstPtr(0);
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
  for (hal_size_t i = 0; i < length; ++i)
  {
    outString[i] = dnaIt.getChar();
    dnaIt.toRight();
  }
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
  for (hal_size_t i = 0; i < length; ++i)
  {
    dnaIt.setChar(inString[i]);
    dnaIt.toRight();
  }
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
      pair<hal_size_t, HDF5Sequence*>(seq->getStartPosition(), seq));
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
  for (i = sequenceDimensions.begin(); i != sequenceDimensions.end(); ++i)
  {
    // Copy segment into HDF5 array
    HDF5Sequence* seq = new HDF5Sequence(this, &_sequenceArray, 
                                         i - sequenceDimensions.begin());
    seq->set(startPosition, *i);
    // Keep the object pointer in our caches
    _sequencePosCache.insert(
      pair<hal_size_t, HDF5Sequence*>(startPosition, seq));
    _sequenceNameCache.insert(pair<string, HDF5Sequence*>(i->_name, seq));
    startPosition += i->_length;
  }  
}
