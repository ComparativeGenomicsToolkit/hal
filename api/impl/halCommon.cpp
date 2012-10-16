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

void hal::reverseComplement(std::string& s)
{
  if (!s.empty())
  {
    size_t j = s.length() - 1;
    size_t i = 0;
    char buf;
    do
    {
      while (j > 0 && s[j] == '-')
      {
        --j;
      }
      while (i < s.length() - 1 && s[i] == '-')
      {
        ++i;
      }
      
      if (i >= j || s[i] == '-' || s[j] == '-')
      {
        if (i == j && s[i] != '-')
        {
          s[i] = reverseComplement(s[i]);
        }
        break;
      }

      buf = reverseComplement(s[i]);
      s[i] = reverseComplement(s[j]);
      s[j] = buf;

      ++i;
      --j;
    } while (true);
  }
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
