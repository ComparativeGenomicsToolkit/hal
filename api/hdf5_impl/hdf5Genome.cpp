/* Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include <iostream>
#include "H5Cpp.h"
#include "hdf5Genome.h"
#include "hdf5DNA.h"
#include "hdf5TopSegment.h"
#include "hdf5BottomSegment.h"

using namespace hal;
using namespace std;
using namespace H5;

const string HDF5Genome::dnaArrayName = "DNA_ARRAY";
const string HDF5Genome::topArrayName = "TOP_ARRAY";
const string HDF5Genome::bottomArrayName = "BOTTOM_ARRAY";
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
}

void HDF5Genome::reset(hal_size_t totalSequenceLength,
                       hal_size_t numTopSegments,
                       hal_size_t numBottomSegments)
{
  H5::Exception::dontPrint();
  try
  {
    DataSet d = _group.openDataSet(dnaArrayName);
    _group.unlink(dnaArrayName);
  }
  catch (H5::Exception){}
  if (totalSequenceLength > 0)
  {
    _dnaArray.create(&_group, dnaArrayName, HDF5DNA::dataType(), 
                     totalSequenceLength, _dcprops);
  }
  resetTopSegments(numTopSegments);
  resetBottomSegments(numBottomSegments);
  _parentCache = NULL;
  _childCache.clear();
}

void HDF5Genome::resetTopSegments(hal_size_t numTopSegments)
{
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

void HDF5Genome::resetBottomSegments(hal_size_t numBottomSegments)
{
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

const string& HDF5Genome::getName() const
{
  return _name;
}

hal_size_t HDF5Genome::getSequenceLength() const
{
  return _dnaArray.getSize();
}

hal_size_t HDF5Genome::getNumberTopSegments() const
{
  return _topArray.getSize();
}

hal_size_t HDF5Genome::getNumberBottomSegments() const
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
   
MetaData* HDF5Genome::getMetaData()
{
  return _metaData;
}

const MetaData* HDF5Genome::getMetaData() const
{
  return _metaData;
}
  
void HDF5Genome::write()
{
  _dnaArray.write();
  _topArray.write();
  _bottomArray.write();
  _metaData->write();
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
}
