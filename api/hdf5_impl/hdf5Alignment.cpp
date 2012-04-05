/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hdf5Alignment.h"
#include "hdf5MetaData.h"

using namespace hal;
using namespace std;
using namespace H5;

/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /MetaData */
const H5std_string HDF5Alignment::MetaGroupName = "MetaData";

HDF5Alignment::HDF5Alignment() :
  Alignment(),
  _file(NULL),
  _cprops(FileCreatPropList::DEFAULT),
  _aprops(FileAccPropList::DEFAULT),
  _flags(H5F_ACC_RDONLY)
{
  
}

HDF5Alignment::~HDF5Alignment()
{
  delete _file;
}

void HDF5Alignment::createNew(const string& alignmentPath)
{
  _flags = H5F_ACC_TRUNC;
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
}

void HDF5Alignment::open(const string& alignmentPath, bool readOnly)
{
  delete _file;
  int _flags = readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
  
}
   
GenomePtr HDF5Alignment::addGenome(const string& path, 
                             const string* parentPath,
                             const vector<string>& childPaths)
{
  return GenomePtr();
}

void HDF5Alignment::removeGenome(const string& path)
{

}

GenomeConstPtr HDF5Alignment::openConstGenome(const string& path) const
{
  return GenomeConstPtr();
}

GenomePtr HDF5Alignment::openGenome(const string& path)
{
  return GenomePtr();
}

string HDF5Alignment::getParent(const string& path)
{
  return "";
}
   
vector<string> HDF5Alignment::getChildren(const string& path)
{
  return vector<string>();
}

MetaDataPtr HDF5Alignment::getMetaData()
{
  return _metaData;
}

MetaDataConstPtr HDF5Alignment::getMetaData() const
{
  return _metaData;
}
