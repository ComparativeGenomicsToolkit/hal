/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include "hdf5Alignment.h"
#include "hdf5MetaData.h"
#include "hdf5Genome.h"
#include "sonLibTree.h"


using namespace hal;
using namespace std;
using namespace H5;

/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /Meta */
const H5std_string HDF5Alignment::MetaGroupName = "Meta";
const H5std_string HDF5Alignment::TreeGroupName = "Phylogeny";
const H5std_string HDF5Alignment::GenomesGroupName = "Genomes";

HDF5Alignment::HDF5Alignment() :
  _file(NULL),
  _cprops(FileCreatPropList::DEFAULT),
  _aprops(FileAccPropList::DEFAULT),
  _dcprops(DSetCreatPropList::DEFAULT),
  _flags(H5F_ACC_RDONLY),
  _metaData(NULL),
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
  _metaData(NULL),
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
  delete _metaData;
  _metaData = new HDF5MetaData(_file, MetaGroupName);
  _tree = NULL;
  _dirty = true;
}

// todo: properly handle readonly
void HDF5Alignment::open(const string& alignmentPath, bool readOnly)
{
  close();
  delete _file;
  int _flags = readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
  _file = new H5File(alignmentPath.c_str(),  _flags, _cprops, _aprops);
  delete _metaData;
  _metaData = new HDF5MetaData(_file, MetaGroupName);
  loadTree();
}
   
void HDF5Alignment::close()
{
  if (_file != NULL)
  {
    writeTree();
    if (_tree != NULL)
    {
      stTree_destruct(_tree);
      _tree = NULL;
    }
    // todo: make sure there's no memory leak with metadata 
    // smart pointer should prevent
    if (_metaData != NULL)
    {
      _metaData->write();
      delete _metaData;
      _metaData = NULL;
    }
  
    map<string, HDF5Genome*>::iterator mapIt;
    for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt)
    {
      HDF5Genome* genome = mapIt->second;
      genome->write();
      delete genome;
    }
    _openGenomes.clear();
    _file->close();
    delete _file;
    _file = NULL;
  }
  else
  {
    assert(_tree == NULL);
    assert(_openGenomes.empty() == true);
  }
}

Genome*  HDF5Alignment::addLeafGenome(const string& name,
                                      const string& parentName,
                                      double branchLength)
{
  if (name.empty() == true || parentName.empty())
  {
    throw hal_exception("name can't be empty");
  }
  map<string, stTree*>::iterator findIt = _nodeMap.find(name);
  if (findIt != _nodeMap.end())
  {
    throw hal_exception(string("node ") + name + " already exists");
  }
  findIt = _nodeMap.find(parentName);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("parent ") + parentName + " not found in tree");
  }
  stTree* parent = findIt->second;
  stTree* node = stTree_construct();
  stTree_setLabel(node, name.c_str());
  stTree_setParent(node, parent);
  stTree_setBranchLength(node, branchLength);
  _nodeMap.insert(pair<string, stTree*>(name, node));

  HDF5Genome* genome = new HDF5Genome(name, this, _file, _dcprops);
  _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  return genome;
}

Genome* HDF5Alignment::addRootGenome(const string& name,
                                        double branchLength)
{
  if (name.empty() == true)
  {
    throw hal_exception("name can't be empty");
  }
  map<string, stTree*>::iterator findIt = _nodeMap.find(name);
  if (findIt != _nodeMap.end())
  {
    throw hal_exception(string("node ") + name + " already exists");
  }
  stTree* node = stTree_construct();
  stTree_setLabel(node, name.c_str());
  if (_tree != NULL)
  {
    stTree_setParent(_tree, node);
    stTree_setBranchLength(_tree, branchLength);
  }
  _tree = node;
  _nodeMap.insert(pair<string, stTree*>(name, node));

  HDF5Genome* genome = new HDF5Genome(name, this, _file, _dcprops);
  _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  return genome;
}


void HDF5Alignment::removeGenome(const string& name)
{
  
}

const Genome* HDF5Alignment::openGenome(const string& name) const
{
  map<string, HDF5Genome*>::iterator mapit = _openGenomes.find(name);
  if (mapit != _openGenomes.end())
  {
    return mapit->second;
  }
  HDF5Genome* genome = new HDF5Genome(name, const_cast<HDF5Alignment*>(this), 
                                      _file, _dcprops);
  genome->read();
  _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  return genome;
}

Genome* HDF5Alignment::openGenome(const string& name)
{
  map<string, HDF5Genome*>::iterator mapit = _openGenomes.find(name);
  if (mapit != _openGenomes.end())
  {
    return mapit->second;
  }
  HDF5Genome* genome = new HDF5Genome(name, this, _file, _dcprops);
  genome->read();
  _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  return genome;
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
  map<string, stTree*>::iterator findIt = _nodeMap.find(name);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("node not found: ") + name);
  }
  stTree* node = findIt->second;
  stTree* parent = stTree_getParent(node);
  if (parent == NULL)
  {
    return "";
  }
  return stTree_getLabel(parent);
}

double HDF5Alignment::getBranchLength(const string& parentName,
                                      const string& childName) const
{
  map<string, stTree*>::iterator findIt = _nodeMap.find(childName);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("node ") + childName + " not found");
  }
  stTree* node = findIt->second;
  stTree* parent = stTree_getParent(node);
  if (parent == NULL || parentName != stTree_getLabel(parent))
  {
    throw hal_exception(string("edge ") + parentName + "--" + childName + 
                        " not found");
  }
  return stTree_getBranchLength(node);
}
   
vector<string> HDF5Alignment::getChildNames(const string& name) const
{
  map<string, stTree*>::iterator findIt = _nodeMap.find(name);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("node ") + name + " not found");
  }
  stTree* node = findIt->second;
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
    assert(_nodeMap.empty() == true);
    return 0;
  }
  else
  {
    return _nodeMap.size();
  }
}

MetaData* HDF5Alignment::getMetaData()
{
  return _metaData;
}

const MetaData* HDF5Alignment::getMetaData() const
{
  return _metaData;
}

string HDF5Alignment::getNewickTree() const
{
  if (_tree == NULL)
  {
    return "";
  }
  else
  {
    char* treeString = stTree_getNewickTreeString(_tree);
    string returnString(treeString);
    free(treeString);
    return returnString;
  }
}

void HDF5Alignment::writeTree()
{
  if (_dirty == false)
     return;

  char* treeString = NULL;
  if (_tree != NULL)
  {
    treeString = stTree_getNewickTreeString(_tree);
  }
  else
  {
    treeString = (char*)malloc(sizeof(char));
    treeString[0] = '\0';
  }
  assert (_file != NULL);
  HDF5MetaData treeMeta(_file, TreeGroupName);
  treeMeta.set(TreeGroupName, treeString);
  free(treeString);
}

static void addNodeToMap(stTree* node, map<string, stTree*>& nodeMap)
{
  const char* label = stTree_getLabel(node);
  assert(label != NULL);
  string name(label);
  assert(nodeMap.find(name) == nodeMap.end());
  nodeMap.insert(pair<string, stTree*>(name, node));
  int32_t numChildren = stTree_getChildNumber(node);
  for (int32_t i = 0; i < numChildren; ++i)
  {
    addNodeToMap(stTree_getChild(node, i), nodeMap);
  }
}

void HDF5Alignment::loadTree()
{
  _nodeMap.clear();
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
    addNodeToMap(_tree, _nodeMap);
  }
}

