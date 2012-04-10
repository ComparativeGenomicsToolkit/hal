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
  map<string, stTree*>::iterator findIt = _nodeMap.find(name);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("node ") + name + " not found");
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
                                      const string& childName)
{
  map<string, stTree*>::iterator findIt = _nodeMap.find(parentName);
  if (findIt == _nodeMap.end())
  {
    throw hal_exception(string("node ") + parentName + " not found");
  }
  stTree* node = findIt->second;
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

static void addNodeToMap(stTree* node, map<string, stTree*> nodeMap)
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

void HDF5Alignment::addGenomeToTree(const string& name,
                                    const pair<string, double>& parent,
                                    const vector<pair<string, double> >&
                                    children)
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

  // construct the new node to add
  stTree* node = stTree_construct();
  stTree_setLabel(node, name.c_str());
  _nodeMap.insert(pair<string, stTree*>(name, node));

  // attach as new root
  if (parent.first.empty() == true)
  {
    stTree_setParent(_tree, node);
    _tree = node;
  }
  
  //attach as new child
  else
  {
    findIt = _nodeMap.find(parent.first);
    if (findIt == _nodeMap.end())
    {
      throw hal_exception(string("parent ") + parent.first + " not found");
    }
    stTree* parentNode = findIt->second;
    stTree_setParent(node, parentNode);
    stTree_setBranchLength(node, parent.second);
  }

  // add in all the children
  vector<pair<string, double> >::const_iterator childIt;
  for (childIt = children.begin(); childIt != children.end(); ++childIt)
  {
    if (childIt->first.empty() == true)
    {
      throw hal_exception("name can't be empty");
    }
    // child already in the tree.  all we do is make sure it is
    // indeed a child of node, and if so set the branch length
    if (_nodeMap.find(childIt->first) != _nodeMap.end())
    {
      stTree* existingChild = stTree_findChild(node, childIt->first.c_str());
      if (existingChild == NULL)
      {
         throw hal_exception("Attempting to add child that violates tree");
      }
      stTree_setBranchLength(existingChild, childIt->second);
    }
    // child not in the tree.  create and add
    else
    {
      stTree* child = stTree_construct();
      stTree_setLabel(child, name.c_str());
      stTree_setParent(child, node);
      stTree_setBranchLength(child, childIt->second);
    }
  }
}
