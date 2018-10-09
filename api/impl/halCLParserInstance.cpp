/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halCLParserInstance.h"
#include "hdf5CLParser.h"
#include "mmapCLParser.h"

using namespace std;
using namespace hal;

CLParserPtr hal::halCLParserInstance(bool createOptions)
{
    CLParserPtr clParser(new CLParser());
    HDF5CLParser::defineOptions(clParser, createOptions);
    MmapCLParser::defineOptions(clParser, createOptions);
    return clParser;
}
