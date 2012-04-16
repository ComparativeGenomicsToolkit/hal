/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALCOMMON_H
#define _HALCOMMON_H

#include <map>
#include <string>
#include <vector>

namespace hal {

std::vector<std::string> chopString(const std::string& inString,
                                    const std::string& separator);

inline hal_dna_t reverseComplement(hal_dna_t c)
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

}
#endif

