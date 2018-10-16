/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halCLParserInstance.h"
#include "hdf5Alignment.h"
#include "mmapAlignment.h"

using namespace std;
using namespace hal;

CLParserPtr hal::halCLParserInstance(unsigned mode)
{
    CLParserPtr parser(new CLParser());
    HDF5Alignment::defineOptions(parser, mode);
    MMapAlignment::defineOptions(parser, mode);
#ifdef ENABLE_UDC
  // this can be used by multiple storage formats
    if ((mode & (CREATE_ACCESS | WRITE_ACCESS)) == 0 ) {
      parser->addOption("udcCacheDir", "udc cache path for *input* hal file(s).",
                        "\"\"");
  }
#endif
    return parser;
}
