/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <map>
#include <cassert>
#include <sys/stat.h>
#include "halCommon.h"
#include "halAlignment.h"
#include "halGenome.h"

using namespace std;
using namespace hal;

const hal_index_t hal::NULL_INDEX = static_cast<hal_index_t>(-1);

const char hal::to_upper_map[128] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
    22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
    41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78,
    79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 65,
    66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
    85, 86, 87, 88, 89, 90, 123, 124, 125, 126, 127};

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

void hal::reverseGaps(string& s)
{
  bool hasGap = false; // just to skip unnecessary work
  string flipped(s.length(), '.');
  for (unsigned i = 0; i < s.length(); ++i)
  {
    if (s[i] == '-')
    {
      hasGap = true;
      flipped[s.length() - 1 - i] = '-';
    }
  }
  if (hasGap == true)
  {
    int j = 0;
    for (unsigned i = 0; i < s.length(); ++i)
    {
      if (s[i] != '-')
      {
        for (; flipped[j] == '-'; ++j);
        flipped.at(j++) = s[i];
      }
    }
    swap(s, flipped);
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

vector<const Genome *> hal::getLeafGenomes(const Alignment* alignment)
{
    // Load in all leaves from alignment
    vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
    vector<const Genome *> leafGenomes;
    for (hal_size_t i = 0; i < leafNames.size(); i++) {
        const Genome *genome = alignment->openGenome(leafNames[i]);
        assert(genome != NULL);
        leafGenomes.push_back(genome);
    }
    return leafGenomes;
}

/* is file a URL that requires UDC? */
bool hal::isUrl(const std::string alignmentPath) {
    return (alignmentPath.find("http:") == 0) or (alignmentPath.find("https:") == 0)
        or (alignmentPath.find("ftp:") == 0);
}

/* get the file size from the OS */
size_t hal::getFileStatSize(int fd) {
    struct stat fileStat;
    if (::fstat(fd, &fileStat) < 0) {
        throw hal_errno_exception("stat failed", errno);
    }
    return fileStat.st_size;
}

/* map of character to encoding for both upper and lower case */
const uint8_t hal::dnaPackMap[256]  = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 11, 4, 4, 4, 10, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

/* map of 4 bit encoding to character */
const char hal::dnaUnpackMap[16] = {'a', 't', 'g', 'c', 'n', '\x00', '\x00', '\x00', 'A', 'T', 'G', 'C', 'N', '\x00', '\x00', '\x00'};


