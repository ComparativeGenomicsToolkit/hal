/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <sstream>
#include <map>
#include "halDefs.h"
#include "hal.h"

using namespace std;
using namespace hal;

const hal_index_t hal::NULL_INDEX = static_cast<hal_index_t>(-1);

/** C++ version of strtok */
vector<string> hal::chopString(const string& inString,
                               const string& separator)
{
  vector<string> outVector;
  string::size_type start = 0;
  string::size_type end = 0;

  while ((end = inString.find(separator, start)) != string::npos)
  {
    // todo: filter out empty strings
    outVector.push_back (inString.substr (start, end-start));
    start = end + separator.size();
  }
  if (start < inString.length())
  {
    outVector.push_back(inString.substr(start, string::npos));
  }
  return outVector;
}

static size_t lcaRecursive(const Genome* genome,
                           const set<const Genome*>& inputSet,
                           map<const Genome*, size_t>& table)
{
  size_t score = 0;
  if (inputSet.find(genome) != inputSet.end())
  {
    score = 1;
  }
  for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
  {
    score += lcaRecursive(genome->getChild(i), inputSet, table);
  }
  table.insert(pair<const Genome*, size_t>(genome, score));
  return score;
}

const Genome* hal::getLowestCommonAncestor(const set<const Genome*>& inputSet)
{
  if (inputSet.empty())
     return NULL;
  // dumb algorithm but it's friday and i'm too tired to think
  const Genome* root = *inputSet.begin();
  while (root->getParent() != NULL)
  {
     root = root->getParent();
  }
  map<const Genome*, size_t> table;
  lcaRecursive(root, inputSet, table);
  const Genome* lca = root;
  bool found = false;
  while (!found)
  {
    found = true;
    hal_size_t score = table.find(lca)->second;
    for (hal_size_t i = 0; found && i < lca->getNumChildren(); ++i)
    {
      if (table.find(lca->getChild(i))->second == score)
      {
        lca = lca->getChild(i);
        found = false;
      }
    }
  }
  return lca;
}

static bool spanningRecursive(const Genome* genome, set<const Genome*>& outputSet, 
                              bool below = false)
{
  bool above = false;
  if (outputSet.find(genome) != outputSet.end())
  {
    below = true;
    above = true;
  }
  hal_size_t numChildren = genome->getNumChildren();
  for (hal_size_t i = 0; i < numChildren; ++i)
  {
    bool childAbove = spanningRecursive(genome->getChild(i), outputSet, below);
    above = above || childAbove;
  }

  if (above && below)
  {
    outputSet.insert(genome);
  }
  return above;  
}

void hal::getGenomesInSpanningTree(const set<const Genome*>& inputSet,
                              set<const Genome*>& outputSet)
{
  const Genome* lca = getLowestCommonAncestor(inputSet);
  if (lca== NULL)
     return;
  outputSet = inputSet;
  outputSet.insert(lca);
  spanningRecursive(lca, outputSet);
}


void hal::getGenomesInSubTree(const Genome* root, 
                              set<const Genome*>& outputSet)
{
  outputSet.insert(root);
  hal_size_t numChildren = root->getNumChildren();
  for (hal_size_t i = 0; i < numChildren; ++i)
  {
    getGenomesInSubTree(root->getChild(i), outputSet);
  }
}

bool PositionCache::insert(hal_index_t pos)
{
  IntervalSet::iterator i = _set.lower_bound(pos);
  IntervalSet::iterator j;
  if (i != _set.end() && i->second <= pos)
  {
    assert(i->first >= pos);
    return false;
  }

  //merge to beginning of existing interval
  if (i != _set.end() && i->second == pos + 1)
  {
    --i->second;
  }
  else
  {
    // create new unit interval
    i = _set.insert(pair<hal_index_t, hal_index_t>(pos, pos)).first;
  }
  assert(i->second <= i->first);
  // merge abutting left interval
  if (i != _set.begin())
  {
    j = i;
    --j;
    if (j->first == i->second - 1)
    {
      i->second = j->second;
      assert(i->second <= i->first);
      _set.erase(j);
    }    
  }

  // merge abutting right interval
  j = i;
  ++j;
  if ( j != _set.end() && j->second == i->first + 1)
  {
     j->second = i->second;
     assert(j->second <= j->first);
     _set.erase(i);
  }

  ++_size;
  assert(find(pos) == true);
  return true;
}

bool PositionCache::find(hal_index_t pos) const
{
  IntervalSet::const_iterator i = _set.lower_bound(pos);
  if (i != _set.end() && i->second <= pos)
  {
    return true;
  }
  return false;
}

void PositionCache::clear()
{
  _set.clear();
  _size = 0;
}

// for debugging
bool PositionCache::check() const
{
  hal_size_t size = 0;
  for (IntervalSet::const_iterator i = _set.begin(); i != _set.end(); ++i)
  {
    size += (i->first + 1) - i->second;
    IntervalSet::const_iterator j = i;
    ++j;
    if (j != _set.end())
    { 
      // test overlap
      if (j->second <= i->first || i->second >= j->first)
      {
        return false;
      }
      // test merge
      if (j->second == i->first + 1)
      {
        return false;
      }
    }
  }
  return size == _size;
}

