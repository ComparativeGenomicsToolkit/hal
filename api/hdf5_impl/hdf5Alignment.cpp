/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <deque>
#include <sstream>
#include <algorithm>
#include "halCommon.h"
#include "hdf5Alignment.h"
#include "hdf5MetaData.h"
#include "hdf5Genome.h"
#include "hdf5CLParser.h"
extern "C" {
#include "sonLibTree.h"
}

using namespace hal;
using namespace std;
using namespace H5;

#ifdef ENABLE_UDC
#include "hdf5UDCFuseDriver.h"
static const hid_t UDC_FUSE_DRIVER_ID =  H5FD_udc_fuse_init();
#endif


/** default group name for MetaData attributes, will be a subgroup
 * of the top level of the file, ie /Meta */
const H5std_string HDF5Alignment::MetaGroupName = "Meta";
const H5std_string HDF5Alignment::TreeGroupName = "Phylogeny";
const H5std_string HDF5Alignment::GenomesGroupName = "Genomes";
const H5std_string HDF5Alignment::VersionGroupName = "Verison";

HDF5Alignment::HDF5Alignment() :
  _file(NULL),
  _flags(H5F_ACC_RDONLY),
  _metaData(NULL),
  _tree(NULL),
  _dirty(false)
{
  // set defaults from the command-line parser
  HDF5CLParser defaultOptions(true);  
  defaultOptions.applyToDCProps(_dcprops);
  defaultOptions.applyToAProps(_aprops);
}

HDF5Alignment::HDF5Alignment(const H5::FileCreatPropList& fileCreateProps,
                             const H5::FileAccPropList& fileAccessProps,
                             const H5::DSetCreatPropList& datasetCreateProps) :
  _file(NULL),
  _flags(H5F_ACC_RDONLY),
  _metaData(NULL),
  _tree(NULL),
  _dirty(false)
{
  _cprops.copy(fileCreateProps);
  _aprops.copy(fileAccessProps);
  _dcprops.copy(datasetCreateProps);
}

HDF5Alignment::~HDF5Alignment()
{
  close();
}

void HDF5Alignment::createNew(const string& alignmentPath)
{
  close();
  _flags = H5F_ACC_TRUNC;
  if (!ofstream(alignmentPath.c_str()))
  {
    throw hal_exception("Unable to open " + alignmentPath);
  }
  setFileDriverFromPath(alignmentPath);
  _file = new H5File(alignmentPath.c_str(), _flags, _cprops, _aprops);
  _file->createGroup(MetaGroupName);
  _file->createGroup(TreeGroupName);
  _file->createGroup(GenomesGroupName);
  _file->createGroup(VersionGroupName);
  delete _metaData;
  _metaData = new HDF5MetaData(_file, MetaGroupName);
  _tree = NULL;
  _dirty = true;
  writeVersion();
}

// todo: properly handle readonly
void HDF5Alignment::open(const string& alignmentPath, bool readOnly)
{
  close();
  delete _file;
  int _flags = readOnly ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
  if (!ifstream(alignmentPath.c_str()))
  {
    throw hal_exception("Unable to open " + alignmentPath);
  }
  setFileDriverFromPath(alignmentPath);
  _file = new H5File(alignmentPath.c_str(),  _flags, _cprops, _aprops);
  if (!compatibleWithVersion(getVersion()))
  {
    stringstream ss;
    ss << "HAL API v" << HAL_VERSION << " incompatible with format v" 
       << getVersion() << " HAL file.";
    throw hal_exception(ss.str());
  }
  delete _metaData;
  _metaData = new HDF5MetaData(_file, MetaGroupName);
  loadTree();
}

// todo: properly handle readonly
void HDF5Alignment::open(const string& alignmentPath) const
{
  const_cast<HDF5Alignment*>(this)->open(alignmentPath, true);
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
    writeVersion();
    map<string, HDF5Genome*>::iterator mapIt;
    for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt)
    {
      HDF5Genome* genome = mapIt->second;
      genome->write();
      delete genome;
    }
    _openGenomes.clear();
    _file->flush(H5F_SCOPE_LOCAL);
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

// same as above but don't write anything to disk.
void HDF5Alignment::close() const
{
  if (_file != NULL)
  {
    if (_tree != NULL)
    {
      stTree_destruct(const_cast<HDF5Alignment*>(this)->_tree);
       const_cast<HDF5Alignment*>(this)->_tree = NULL;
    }
    // todo: make sure there's no memory leak with metadata 
    // smart pointer should prevent
    if (_metaData != NULL)
    {
      delete  const_cast<HDF5Alignment*>(this)->_metaData;
       const_cast<HDF5Alignment*>(this)->_metaData = NULL;
    }

    map<string, HDF5Genome*>::iterator mapIt;
    for (mapIt = _openGenomes.begin(); mapIt != _openGenomes.end(); ++mapIt)
    {
      HDF5Genome* genome = mapIt->second;
      delete genome;
    }
    _openGenomes.clear();
     const_cast<HDF5Alignment*>(this)->_file->close();
     delete const_cast<HDF5Alignment*>(this)->_file;
     const_cast<HDF5Alignment*>(this)->_file = NULL;
  }
  else
  {
    assert(_tree == NULL);
    assert(_openGenomes.empty() == true);
  }
}

void HDF5Alignment::setOptionsFromParser(CLParserConstPtr parser) const
{
  const HDF5CLParser* hdf5Parser = 
     dynamic_cast<const HDF5CLParser*>(parser.get());
  if (hdf5Parser == NULL)
  {
    throw hal_exception("setOptionsFromParser expected instance of "
                        "hdf5CLParser");
  }
  hdf5Parser->applyToDCProps(_dcprops);
  hdf5Parser->applyToAProps(_aprops);
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
  _dirty = true;
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
  _dirty = true;
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
  HDF5Genome* genome = NULL;
  if (_nodeMap.find(name) != _nodeMap.end())
  {
    genome = new HDF5Genome(name, const_cast<HDF5Alignment*>(this), 
                                        _file, _dcprops);
    genome->read();
    _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  }
  return genome;
}

Genome* HDF5Alignment::openGenome(const string& name)
{
  map<string, HDF5Genome*>::iterator mapit = _openGenomes.find(name);
  if (mapit != _openGenomes.end())
  {
    return mapit->second;
  }
  HDF5Genome* genome = NULL;
  if (_nodeMap.find(name) != _nodeMap.end())
  {
    genome = new HDF5Genome(name, this, _file, _dcprops);
    genome->read();
    _openGenomes.insert(pair<string, HDF5Genome*>(name, genome));
  }
  return genome;
}

void HDF5Alignment::closeGenome(const Genome* genome) const
{
  string name = genome->getName();
  map<string, HDF5Genome*>::iterator mapIt = _openGenomes.find(name);
  if (mapIt == _openGenomes.end())
  {
    throw hal_exception("Attempt to close non-open genome.  " 
                        "Should not even be possible");
  }
  mapIt->second->write();
  delete mapIt->second;
  _openGenomes.erase(mapIt);
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
  vector<string> leaves;
  vector<string> children;
  deque<string> bfQueue;
  bfQueue.push_front(name);
  while (bfQueue.empty() == false)
  {
    string& current = bfQueue.back();
    children = getChildNames(current);
    if (children.empty() == true && current != name)
    {
      leaves.push_back(current);
    }
    for (size_t i = 0; i < children.size(); ++i)
    {
      bfQueue.push_front(children[i]);
    }
    bfQueue.pop_back();
  }
  return leaves;
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

string HDF5Alignment::getVersion() const
{
  try
  {
    H5::Exception::dontPrint();
    _file->openGroup(VersionGroupName);  
    HDF5MetaData versionMeta(_file, VersionGroupName);
    if (versionMeta.has(VersionGroupName) == false)
    {
      throw Exception();
    }
    return versionMeta.get(VersionGroupName);
  }
  catch (Exception& e)
  {
    // since there was no version tag at the beginning,
    // we retroactively consider it 0.0
    return "0.0";
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

void HDF5Alignment::writeVersion()
{
  if (_dirty == false)
     return;
  
  assert(_file != NULL);
  HDF5MetaData versionMeta(_file, VersionGroupName);
  stringstream ss;
  ss << HAL_VERSION;
  versionMeta.set(VersionGroupName, ss.str());
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

void HDF5Alignment::setFileDriverFromPath(const string& path)
{
#ifdef ENABLE_UDC
  string pathCpy(path);
  transform(pathCpy.begin(), pathCpy.end(), pathCpy.begin(), ::tolower);
  if (pathCpy.find("http:") == 0 || 
      pathCpy.find("https:") == 0 ||
      pathCpy.find("ftp:") == 0)
  {
    _aprops.setDriver(UDC_FUSE_DRIVER_ID, NULL);
  }
#endif
}
