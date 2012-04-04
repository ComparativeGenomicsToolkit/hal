/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <cassert>
#include "hdf5MetaData.h"
#include "halCommon.h"

using namespace hal;
using namespace H5;
using namespace std;

HDF5MetaData::HDF5MetaData() :
  _parent(NULL)
{
}

HDF5MetaData::HDF5MetaData(CommonFG* parent, const string& name)
{
  open(parent, name);
}

HDF5MetaData::~HDF5MetaData()
{
  write();
}

void HDF5MetaData::set(const string& key, const string& value)
{
  if (has(key) == true)
  {
    _map[key] = value;
  }
  else
  {
    _map.insert(pair<string, string>(key, value));
  }
  _dirty = true;
}

const string& HDF5MetaData::get(const string& key) const
{
  assert (has(key) == true);
  return _map.find(key)->second;
}

bool HDF5MetaData::has(const string& key) const
{
  return _map.find(key) != _map.end();
}

const map<string, string>& HDF5MetaData::getMap() const
{
  return _map;
}

static void attr_operator(H5Object& loc/*in*/,
                          const H5std_string attr_name/*in*/,
                          void *operator_data/*in,out*/)
{
  map<string, string>* attMap = 
     static_cast<map<string, string>*>(operator_data);
  StrType vlsType(0, H5T_VARIABLE);
  Attribute attr = loc.openAttribute(attr_name);
  H5std_string strg;
  attr.read(vlsType, strg);
  attMap->insert(pair<string, string>(attr_name, strg));
}

void HDF5MetaData::open(CommonFG* parent, const string& name)
{
  assert(parent != NULL);
  _map.clear();
  _parent = parent;
  _name = name;
  _dirty = false;

  try
  {
    _group = parent->openGroup(name);
  }
  catch (Exception& e)
  {
    _group = parent->createGroup(name);
  }

  if (_group.getNumAttrs() > 0)
  {
    _group.iterateAttrs(attr_operator, NULL, &_map);
  }
  assert(_map.size() == (size_t)_group.getNumAttrs());
}

void HDF5MetaData::write()
{
  if (!_dirty)
     return;

  for (map<string, string>::iterator i = _map.begin();
       i != _map.end(); 
       ++i)
  {
    StrType vlsType(0, H5T_VARIABLE);
    DataSpace attSpace(H5S_SCALAR);
    Attribute attr;
    // we delve into C api for this test
    if (H5Aexists(_group.getId(), i->first.c_str()) == true)
    {
      _group.removeAttr(_name);
    }    
    attr = _group.createAttribute(i->first, vlsType, attSpace);
    attr.write(vlsType, i->second);
  }
  _dirty = false;
}
