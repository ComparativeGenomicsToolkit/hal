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

// we now work with names instead of Genome*s to avoid expensive openGenome
// function
static size_t lcaRecursive(const Alignment* alignment,
                           const string& genome,
                           const set<string>& inputSet,
                           map<string, size_t>& table)
{
  size_t score = 0;
  if (inputSet.find(genome) != inputSet.end())
  {
    score = 1;
  }
  vector<string> childs = alignment->getChildNames(genome);
  for (hal_size_t i = 0; i < childs.size(); ++i)
  {
    score += lcaRecursive(alignment, childs[i], inputSet, table);
  }
  table.insert(pair<string, size_t>(genome, score));
  return score;
}

const Genome* hal::getLowestCommonAncestor(const set<const Genome*>& inputSet)
{
  if (inputSet.empty())
     return NULL;
  
  const Alignment* alignment = (*inputSet.begin())->getAlignment();
  set<string> inputNames;
  for (set<const Genome*>::iterator i = inputSet.begin(); i != inputSet.end();
       ++i)
  {
    inputNames.insert((*i)->getName());
  }

  // dumb algorithm but it's friday and i'm too tired to think
  string root = alignment->getRootName();
  map<string, size_t> table;
  lcaRecursive(alignment, root, inputNames, table);
  string lca = root;
  bool found = false;
  while (!found)
  {
    found = true;
    hal_size_t score = table.find(lca)->second;
    vector<string> childs = alignment->getChildNames(lca);

    for (hal_size_t i = 0; found && i < childs.size(); ++i)
    {
      if (table.find(childs[i])->second == score)
      {
        lca = childs[i];
        found = false;
      }
    }
  }
  return alignment->openGenome(lca);
}

// we now work with names instead of Genome*s to avoid expensive openGenome
// function
static bool spanningRecursive(const Alignment* alignment,
                              const string& genome, 
                              set<const Genome*>& outputSet, 
                              set<string>& outputNames,
                              bool below = false)
{
  bool above = false;
  if (outputNames.find(genome) != outputNames.end())
  {
    below = true;
    above = true;
  }
  vector<string> childs = alignment->getChildNames(genome);
  for (hal_size_t i = 0; i < childs.size(); ++i)
  {
    bool childAbove = spanningRecursive(alignment, childs[i], outputSet, 
                                        outputNames, below);
    above = above || childAbove;
  }

  if (above && below)
  {
    outputSet.insert(alignment->openGenome(genome));
    outputNames.insert(genome);
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
  set<string> outputNames;
  for (set<const Genome*>::iterator i = outputSet.begin(); i != outputSet.end();
       ++i)
  {
    outputNames.insert((*i)->getName());
  }
  spanningRecursive(lca->getAlignment(), lca->getName(), outputSet, outputNames);
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
