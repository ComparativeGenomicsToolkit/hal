/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "hdf5Alignment.h"
#include "hdf5MetaData.h"
#include "sonLibTree.h"


using namespace hal;
using namespace std;
using namespace H5;

/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /MetaData */
const H5std_string HDF5Alignment::MetaGroupName = "__MetaData";
const H5std_string HDF5Alignment::TreeGroupName = "__PhylogenyData";

HDF5Alignment::HDF5Alignment() :
  _file(NULL),
  _cprops(FileCreatPropList::DEFAULT),
  _aprops(FileAccPropList::DEFAULT),
  _dcprops(DSetCreatPropList::DEFAULT),
  _flags(H5F_ACC_RDONLY),
  _tree(NULL)
{
  
}


HDF5Alignment::HDF5Alignment(const H5::FileCreatPropList& fileCreateProps,
                             const H5::FileAccPropList& fileAccessProps,
                             const H5::DSetCreatPropList& datasetCreateProps) :
  _file(NULL),
  _cprops(fileCreateProps),
  _aprops(fileAccessProps),
  _dcprops(datasetCreateProps),
  _flags(H5F_ACC_RDONLY),
  _tree(NULL)
{

}

HDF5Alignment::~HDF5Alignment()
{
  close();
}

void HDF5Alignment::createNew(const string& alignmentPath)
{
  _flags = H5F_ACC_TRUNC;
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
  _tree = stTree_construct();
}

void HDF5Alignment::open(const string& alignmentPath, bool readOnly)
{
  delete _file;
  int _flags = readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
  
}
   
void HDF5Alignment::close()
{
  if (_file != NULL)
  {
     _file->close();
     delete _file;
     _file = NULL;
     assert(_tree != NULL);
  }
  if (_tree != NULL)
  {
    writeTree();
    stTree_destruct(_tree);
    _tree = NULL;
  }
}

GenomePtr HDF5Alignment::addGenome(const string& name,
                                   const string& path,
                                   const string& parentName,
                                   const vector<string>& childNames)
{
  return GenomePtr();
}

void HDF5Alignment::removeGenome(const string& name)
{

}

GenomeConstPtr HDF5Alignment::openConstGenome(const string& name) const
{
  return GenomeConstPtr();
}

GenomePtr HDF5Alignment::openGenome(const string& name)
{
  return GenomePtr();
}

string HDF5Alignment::getRootName() const
{
  return "";
}

string HDF5Alignment::getParentName(const string& name) const
{
  return "";
}
   
vector<string> HDF5Alignment::getChildNames(const string& name) const
{
  return vector<string>();
}

vector<string> HDF5Alignment::getLeafNamesBelow(const string& name) const
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

const string& HDF5Alignment::getPath(const string& name) const
{
  return "";
}
