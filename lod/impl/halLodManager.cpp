/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <deque>
#include <sstream>
#include <limits>
#include <fstream>
#include <algorithm>
#include "halLodManager.h"

#ifdef ENABLE_UDC
extern "C" {
#include "common.h"
#include "udc.h"
// Dont want this macro ruining std::max
#ifdef max
#undef max
#endif
}
#endif

using namespace std;
using namespace hal;

// default to 5 days for now
const unsigned long LodManager::MaxAgeSec = 432000;

// specify upper limit of lods.
// (MUST MANUALLY KEEP CONSISTENT WITH global MaxLodToken variable in 
// hal/lod/halLodInterpolate.py)
const string LodManager::MaxLodToken = "max";

LodManager::LodManager()
{

}

LodManager::~LodManager()
{
  for (AlignmentMap::iterator mapIt = _map.begin(); mapIt != _map.end();
       ++mapIt)
  {
    if (mapIt->second.second.get() != NULL)
    {
      mapIt->second.second->close();
    }
  }
}

void LodManager::loadLODFile(const string& lodPath,
                             CLParserConstPtr options)
{
  _map.clear();
  bool loadAll = false;

#ifdef ENABLE_UDC
  char* cpath = const_cast<char*>(lodPath.c_str());
  if (lodPath.find("http") == 0)
  {
    unsigned long cpathAge = udcCacheAge(cpath, NULL);
    loadAll = cpathAge > MaxAgeSec;
  }
  
  size_t cbufSize = 0;
  char* cbuffer = udcFileReadAll(cpath, NULL, 100000, &cbufSize);
  if (cbuffer == NULL)
  {
    stringstream emes;
    emes << "Error udc-opening " << lodPath;
    throw hal_exception(emes.str());
  }
  string cbufCpy(cbuffer);
  freeMem(cbuffer);
  stringstream ifile(cbufCpy);
#else
  ifstream ifile(lodPath.c_str());
#endif

  if (!ifile.good())
  {
    stringstream emes;
    emes << "Error opening " << lodPath;
    throw hal_exception(emes.str());
  }
  
  string lineBuffer;
  hal_size_t minLen;
  string path;
  hal_size_t lineNum = 1;
  bool foundMax = false;
  _maxLodLowerBound = (hal_size_t)numeric_limits<hal_index_t>::max();
  while (ifile.good())
  {
    getline(ifile, lineBuffer);
    stringstream ss(lineBuffer);
    ss >> minLen >> path;
    if (ifile.bad())
    {
      stringstream emes;
      emes << "Error parsing line " << lineNum << " of " << lodPath;
      throw hal_exception(emes.str());
    }
    if (lineBuffer.length() == 0)
    {
      continue;
    }
    string fullHalPath;
    if (foundMax == true)
    {
      stringstream emes;
      emes << "Error on line " << lineNum << " of " << lodPath 
           << ": Limit token (" << MaxLodToken << ") can only appear on"
           << " final line.";
        throw hal_exception(emes.str());
    }
    if (path == MaxLodToken)
    {
      foundMax = true;
      _maxLodLowerBound = minLen;
      fullHalPath = MaxLodToken;
    }
    else
    {
      fullHalPath = resolvePath(lodPath, path);
    }
    _map.insert(pair<hal_size_t, PathAlign>(
                  minLen, PathAlign(fullHalPath, AlignmentConstPtr())));
    ++lineNum;
  }

  checkMap(lodPath);
  if (loadAll == true)
  {
    preloadAlignments();
  }
}

void LodManager::loadSingeHALFile(const string& halPath,
                                  CLParserConstPtr options)
{
  _map.clear();
  _map.insert(pair<hal_size_t, PathAlign>(
                0, PathAlign(halPath, AlignmentConstPtr())));
  _maxLodLowerBound = (hal_size_t)numeric_limits<hal_index_t>::max();
  checkMap(halPath);
}

AlignmentConstPtr LodManager::getAlignment(hal_size_t queryLength,
                                           bool needDNA)
{
  assert(_map.size() > 0);
  AlignmentMap::iterator mapIt;
  if (needDNA == true)
  {
    mapIt = _map.begin();
  }
  else
  {
    mapIt = _map.upper_bound(queryLength);
    --mapIt;
  }
  assert(mapIt->first <= queryLength);
  AlignmentConstPtr& alignment = mapIt->second.second;
  if (mapIt->first == _maxLodLowerBound)
  {
    stringstream ss;
    ss << "Query length " << queryLength << " above maximum LOD size of "
       << getMaxQueryLength();
    throw hal_exception(ss.str());
  }
  if (alignment.get() == NULL)
  {
    alignment = hdf5AlignmentInstanceReadOnly();
    if (_options.get() != NULL)
    {
      alignment->setOptionsFromParser(_options);
    }
    alignment->open(mapIt->second.first);
    checkAlignment(mapIt->first, mapIt->second.first, alignment);
  }
  assert(mapIt->second.second.get() != NULL);
  return alignment;
}

bool LodManager::isLod0(hal_size_t queryLength) const
{
  assert(_map.size() > 0);
  AlignmentMap::const_iterator mapIt = _map.upper_bound(queryLength);
  --mapIt;
  return mapIt == _map.begin();
}

string LodManager::resolvePath(const string& lodPath,
                               const string& halPath)
{
  assert(lodPath.empty() == false && halPath.empty() == false);
  if (halPath[0] == '/' || halPath.find(":/") != string::npos)
  {
    return halPath;
  }
  size_t sPos = lodPath.find_last_of('/');
  if (sPos == string::npos)
  {
    return halPath;
  }
  return lodPath.substr(0, sPos + 1) + halPath;
}

void LodManager::checkMap(const string& lodPath)
{
  if (_map.size() == 0)
  {
    stringstream ss;
    ss << "No entries were found in " << lodPath;
    throw hal_exception(ss.str());
  }
  AlignmentMap::const_iterator mapIt = _map.begin();
  if (mapIt->first != 0)
  {
    stringstream ss;
    ss << "No alignment with range value 0 found in " << lodPath << ". "
       << "A record of the form \"0 pathToOriginalHALFile\" must be present";
    throw hal_exception(ss.str());
  }
  if (_maxLodLowerBound == 0)
  {
    throw hal_exception("Maximum LOD query length must be > 0");
  }
}

void LodManager::checkAlignment(hal_size_t minQuery,
                                const string& path,
                                AlignmentConstPtr alignment)
{
  if (alignment->getNumGenomes() == 0)
  {
    stringstream ss;
    ss << "No genomes found in base alignment specified in " << path;
    throw hal_exception(ss.str());
  }

#ifndef NDEBUG
  if (minQuery == 0)
  {
    vector<string> leafNames = alignment->getLeafNamesBelow(
      alignment->getRootName());
    string name = !leafNames.empty() ? leafNames[0] : alignment->getRootName();
    const Genome* genome = alignment->openGenome(name);
    
    bool seqFound = genome->containsDNAArray();
    alignment->closeGenome(genome);
    if (seqFound == false)
    {
      stringstream ss;
      ss << "HAL file for highest level of detail (0) in " << path 
         << "must contain DNA sequence information.";
      throw hal_exception(ss.str());
    }
  }
#endif
}

void LodManager::preloadAlignments()
{
  for (AlignmentMap::iterator i = _map.begin(); i != _map.end(); ++i)
  {
    AlignmentConstPtr alignment = getAlignment(i->first, false);
    if (alignment->getNumGenomes() > 0)
    {
      const Genome* root = alignment->openGenome(alignment->getRootName()); 
      (void)root;
//      set<const Genome*> genomeSet;
//      getGenomesInSubTree(root, genomeSet);
     }
  }
}

