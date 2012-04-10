/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include "hdf5Alignment.h"
#include "hdf5MetaData.h"
#include "hdf5Genome.h"
#include "sonLibTree.h"


using namespace hal;
using namespace std;
using namespace H5;

/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /Meta */
const H5std_string HDF5Alignment::MetaGroupName = "/Meta";
const H5std_string HDF5Alignment::TreeGroupName = "/Phylogeny";
const H5std_string HDF5Alignment::GenomesGroupName = "/Genomes";

HDF5Alignment::HDF5Alignment() :
  _file(NULL),
  _cprops(FileCreatPropList::DEFAULT),
  _aprops(FileAccPropList::DEFAULT),
  _dcprops(DSetCreatPropList::DEFAULT),
  _flags(H5F_ACC_RDONLY),
  _tree(NULL),
  _dirty(false)
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
  _tree(NULL),
  _dirty(false)
{

}

HDF5Alignment::~HDF5Alignment()
{
  close();
}

void HDF5Alignment::createNew(const string& alignmentPath)
{
  close();
  _flags = H5F_ACC_TRUNC;
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _file->createGroup(MetaGroupName);
  _file->createGroup(TreeGroupName);
  _file->createGroup(GenomesGroupName);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
  _tree = NULL;
}

// todo: properly handle readonly
void HDF5Alignment::open(const string& alignmentPath, bool readOnly)
{
  close();
  delete _file;
  int _flags = readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
  _file = new H5File(alignmentPath.c_str(),  _flags, _cprops, _aprops);
  _metaData = MetaDataPtr(new HDF5MetaData(_file, MetaGroupName));
  loadTree();
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
  // todo: make sure there's no memory leak with metadata 
  // smart pointer should prevent
  dynamic_cast<HDF5MetaData*>(_metaData.get())->write();
  
  map<string, GenomePtr>::iterator mapIt;
  for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt)
  {
    HDF5Genome* genome = dynamic_cast<HDF5Genome*>(mapIt->second.get());
    genome->write();
  }
  _openGenomes.clear();
}

GenomePtr HDF5Alignment::addGenome(const string& name,
                                   const pair<string, double>& parent,
                                   const vector<pair<string, double> >& 
                                   children)
{
  addGenomeToTree(name, parent, children);
  HDF5Genome* genome = new HDF5Genome(name, this, _file, _dcprops);
  GenomePtr genomePtr(genome);
  _openGenomes.insert(pair<string, GenomePtr>(name, genomePtr));
  return genomePtr;
}

void HDF5Alignment::removeGenome(const string& name)
{
  
}

GenomeConstPtr HDF5Alignment::openConstGenome(const string& name) const
{
  map<string, GenomePtr>::iterator mapit = _openGenomes.find(name);
  if (mapit != _openGenomes.end())
  {
    return mapit->second;
  }
  HDF5Genome* genome = new HDF5Genome(name, const_cast<HDF5Alignment*>(this), 
                                      _file, _dcprops);
  genome->read();
  GenomePtr genomePtr(genome);
  _openGenomes.insert(pair<string, GenomePtr>(name, genomePtr));
  return genomePtr;
}

GenomePtr HDF5Alignment::openGenome(const string& name)
{
  map<string, GenomePtr>::iterator mapit = _openGenomes.find(name);
  if (mapit != _openGenomes.end())
  {
    return mapit->second;
  }
  HDF5Genome* genome = new HDF5Genome(name, this, _file, _dcprops);
  genome->read();
  GenomePtr genomePtr(genome);
  _openGenomes.insert(pair<string, GenomePtr>(name, genomePtr));
  return genomePtr;
}

string HDF5Alignment::getRootName() const
{
  if (_tree == NULL)
  {
    throw hal_exception("Can't get root name of empty tree");
  }
  return stTree_getLabel(_tree);
}

string HDF5Alignment::getParentName(const string& name) const
{
  stTree* node = stTree_findChild(_tree, name.c_str());
  if (node == NULL)
  {
    throw hal_exception(string("node ") + name + " not found");
  }
  stTree* parent = stTree_getParent(node);
  if (parent == NULL)
  {
    return "";
  }
  return stTree_getLabel(parent);
}

double HDF5Alignment::getBranchLength(const string& parentName,
                                      const string& childName)
{
  stTree* node = stTree_findChild(_tree, childName.c_str());
  if (node == NULL)
  {
    throw hal_exception(string("node ") + childName + " not found");
  }
  stTree* parent = stTree_getParent(node);
  if (parent == NULL || parentName != stTree_getLabel(parent))
  {
    throw hal_exception(string("edge ") + childName + "--" + parentName + 
                        " not found");
  }
  return stTree_getBranchLength(node);
}
   
vector<string> HDF5Alignment::getChildNames(const string& name) const
{
  stTree* node = stTree_findChild(_tree, name.c_str());
  if (node == NULL)
  {
    throw hal_exception(string("node ") + name + " not found");
  }
  int32_t numChildren = stTree_getChildNumber(node);
  vector<string> childNames(numChildren);
  for (int32_t i = 0; i < numChildren; ++i)
  {
    childNames[i] = stTree_getLabel(stTree_getChild(node, i));
  }
  return childNames;
}

vector<string> HDF5Alignment::getLeafNamesBelow(const string& name) const
{                                               
  return vector<string>();
}

hal_size_t HDF5Alignment::getNumGenomes() const
{
  if (_tree == NULL)
  {
    return 0;
  }
  else
  {
    return stTree_getNumNodes(_tree);
  }
}

MetaDataPtr HDF5Alignment::getMetaData()
{
  return _metaData;
}

MetaDataConstPtr HDF5Alignment::getMetaData() const
{
  return _metaData;
}

void HDF5Alignment::writeTree()
{
  if (_dirty == false)
     return;

  assert(_tree != NULL);
  
  char* tree = stTree_getNewickTreeString(_tree);
  HDF5MetaData treeMeta(_file, TreeGroupName);
  treeMeta.set(TreeGroupName, tree);
  free(tree);
}

void HDF5Alignment::loadTree()
{
  HDF5MetaData treeMeta(_file, TreeGroupName);
  map<string, string> metaMap = treeMeta.getMap();
  assert(metaMap.size() == 1);
  assert(metaMap.find(TreeGroupName) != metaMap.end());
  const string& treeString = metaMap[TreeGroupName];
  if (_tree != NULL)
  {
    stTree_destruct(_tree);
  }
  if (treeString.empty() == true)
  {
    _tree = stTree_construct();
  }
  else
  {
    _tree = stTree_parseNewickString(const_cast<char*>(treeString.c_str()));
  }
}
