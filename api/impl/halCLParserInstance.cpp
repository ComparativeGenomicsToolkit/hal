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

CLParserPtr hal::halCLParserInstance(unsigned mode)
{
    CLParserPtr parser(new CLParser());
    HDF5CLParser::defineOptions(parser, mode);
    MMapCLParser::defineOptions(parser, mode);
#ifdef ENABLE_UDC
  // this can be define by other storage formats as well
    if ((mode & (CREATE_ACCESS | WRITE_ACCESS)) == 0 ) {
      parser->addOption("udcCacheDir", "udc cache path for *input* hal file(s).",
                        "\"\"");
  }
#endif
    return parser;
}
