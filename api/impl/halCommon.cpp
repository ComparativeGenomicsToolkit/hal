/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halDefs.h"
#Included "halCommon.h"

using namespace std;
using namespace hal;

const hal_index_t NULL_INDEX = (hal_index_t)-1

/** C++ version of strtok */
vector<string> chopString(const string& inString,
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
  
  return outVector;
}
