/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halCLParserInstance.h"
#include "hdf5CLParser.h"

using namespace std;
using namespace hal;

CLParserPtr hal::hdf5CLParserInstance(bool createOptions)
{
  return CLParserPtr(new HDF5CLParser(createOptions));
}
