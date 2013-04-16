/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <deque>
#include <sstream>
#include <fstream>
#include "halLodManager.h"
 
using namespace std;
using namespace hal;

LodManager::LodManager()
{

}

LodManager::~LodManager()
{
  for (AlignmentMap::iterator mapIt = _map.begin(); mapIt != _map.end();
       ++mapIt)
  {
    mapIt->second->close();
  }
}

void LodManager::loadLODFile(const string& lodPath,
                             CLParserConstPtr options)
{
  _map.clear();
  ifstream ifile(lodPath.c_str());
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
    AlignmentConstPtr alignment = hdf5AlignmentInstanceReadOnly();
    if (options.get() != NULL)
    {
      alignment->setOptionsFromParser(options);
    }
    alignment->open(path);
    _map.insert(pair<hal_size_t, AlignmentConstPtr>(minLen, alignment));
    ++lineNum;
  }
  checkMap(lodPath);
}

void LodManager::loadSingeHALFile(const string& halPath,
                                  CLParserConstPtr options)
{
  _map.clear();
  AlignmentConstPtr alignment = hdf5AlignmentInstanceReadOnly();
  if (options.get() != NULL)
  {
    alignment->setOptionsFromParser(options);
  }
  alignment->open(halPath);
  _map.insert(pair<hal_size_t, AlignmentConstPtr>(0, alignment));
  checkMap(halPath);
}

AlignmentConstPtr LodManager::getAlignment(hal_size_t queryLength) const
{
  assert(_map.size() > 0);
  AlignmentMap::const_iterator mapIt = _map.upper_bound(queryLength);
  --mapIt;
  assert(mapIt->first <= queryLength);
  return mapIt->second;
}

void LodManager::checkMap(const string& lodPath) const
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
  AlignmentConstPtr alignment = mapIt->second;
  if (alignment->getNumGenomes() == 0)
  {
    stringstream ss;
    ss << "No genomes found in base alignment specified in " << lodPath;
    throw hal_exception(ss.str());
  }

  bool seqFound = false;
  deque<string> bfQueue;
  bfQueue.push_back(alignment->getRootName());
  while (bfQueue.size() > 0 && !seqFound)
  {
    string name = bfQueue.back();
    bfQueue.pop_back();
    const Genome* genome = alignment->openGenome(name);
    seqFound = genome->containsDNAArray();
    alignment->closeGenome(genome);
    vector<string> children = alignment->getChildNames(name);
    for (size_t i = 0; i < children.size(); ++i)
    {
      bfQueue.push_front(children[i]);
    }
  }
  if (seqFound == false)
  {
    stringstream ss;
    ss << "HAL file for highest level of detail (0) in " << lodPath 
       << "must contain DNA sequence information.";
    throw hal_exception(ss.str());
  }
}
