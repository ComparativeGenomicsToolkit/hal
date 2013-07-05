/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOMMON_H
#define _HALCOMMON_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <locale>
#include <cassert>
#include <sstream>
#include "hal.h"

namespace hal {

inline bool compatibleWithVersion(const std::string& version)
{
  double myVersion, inVersion;
  // assume versions are strings tho we treat as floats for now.
  std::stringstream ss, ss2;
  ss << HAL_VERSION;
  ss >> myVersion;
  ss2 << version;
  ss2 >> inVersion;
  return (int)myVersion == (int)inVersion;
}

/** C++ style strtok-type function.  Can't remember why I wrote it */
std::vector<std::string> chopString(const std::string& inString,
                                    const std::string& separator);

/** Get the DNA reverse complement of a character.
 * If the input is not a nucleotide, then just return it as is
 * (ie no error checking) */
inline char reverseComplement(char c)
{
  switch (c)
  {
  case 'A' : return 'T'; 
  case 'a' : return 't'; 
  case 'C' : return 'G'; 
  case 'c' : return 'g';
  case 'G' : return 'C';
  case 'g' : return 'c';
  case 'T' : return 'A';
  case 't' : return 'a';
  default : break;
  }
  return c;
}

/** Get the reversed complement of a string (in place */
void reverseComplement(std::string& s);

/** Check if a DNA character is a valid base (or n-chracter) */
inline bool isNucleotide(char c)
{
  bool result = false;
  switch (c)
  {
  case 'A' : 
  case 'a' : 
  case 'C' : 
  case 'c' : 
  case 'G' : 
  case 'g' : 
  case 'T' : 
  case 't' : 
  case 'N' :
  case 'n' :
    result = true;
  default : break;
  }
  return result;
}

inline bool isTransition(char c1, char c2)
{
  assert(isNucleotide(c1) && isNucleotide(c2));
  char x = std::toupper((char)c1);
  char y = std::toupper((char)c2);
  switch(x)
  {
  case 'A' : return y == 'G';
  case 'C' : return y == 'T';
  case 'G' : return y == 'A';
  case 'T' : return y == 'C';
  default: break;
  }
  return false;
}

inline bool isSubstitution(char c1, char c2)
{
  return std::toupper(c1) != std::toupper(c2);
}

inline bool isTransversion(char c1, char c2)
{
  char x = std::toupper((char)c1);
  char y = std::toupper((char)c2);
  return (x != y && x != 'N' && y != 'N' && !isTransition(c1, c2));
}

inline bool isMissingData(char c)
{
  return c == 'n' || c == 'N';
}

inline bool isMasked(char c)
{
  return c == std::tolower(c);
}

/** test if 3rd codon position is 4-fold degenerate given first 2 positions */
inline bool isFourfoldDegenerate(char c1, char c2)
{
  char x1 = std::toupper((char)c1);
  char x2 = std::toupper((char)c2);
  if (x2 == 'T' || x2 == 'G')
  {
    return x1 == 'C' || x1 == 'G';
  }
  else if (x2 == 'C')
  {
    return x1 == 'A' || x1 == 'C' || x1 == 'G' || x1 == 'T';
  }
  return false;
}

/** Count the mutations between two DNA strings */
inline hal_size_t hammingDistance(const std::string& s1, const std::string& s2)
{
  assert(s1.length() == s2.length());
  hal_size_t dist = 0;
  for (size_t i = 0; i < s1.length(); ++i)
  {
    if (isSubstitution(s1[i], s2[i]) == true)
    {
      ++dist;
    }
  }
  return dist;
}

const Genome* getLowestCommonAncestor(const std::set<const Genome*>& inputSet);

/* Given a set of genomes (input set) find all genomes in the spanning
 * tree including the inptuts (root should be the root of the alignment) */
void getGenomesInSpanningTree(const std::set<const Genome*>& inputSet,
                              std::set<const Genome*>& outputSet);

/* Given a node (root), return it and all genomes (including internal nodes)
 * below it in the tree */
void getGenomesInSubTree(const Genome* root, 
                         std::set<const Genome*>& outputSet);

}

#endif
