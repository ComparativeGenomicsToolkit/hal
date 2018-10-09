/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "mmapCLParser.h"
#include "mmapFile.h"

using namespace hal;
using namespace std;

void MmapCLParser::defineOptions(CLParserPtr parser,
                                 bool createOptions)
{
  if (createOptions)
  {
      parser->addOption("mmapInitSize", "mmap HAL file initial size", MMAP_DEFAULT_INIT_SIZE);
      parser->addOption("mmapGrowSize", "mmap HAL file size to grow when more memory is needed", MMAP_DEFAULT_INIT_SIZE);
  }
#ifdef ENABLE_UDC
  // this can be define by other storage formats as well
  if (not parser->hasArgument("udcCacheDir")) {
      parser->addOption("udcCacheDir", "udc cache path for *input* hal file(s).",
                        "\"\"");
  }
#endif
}

size_t MmapCLParser::getInitSize(CLParserPtr parser)
{
    return parser->get<size_t>("mmapInitSize");
}

size_t MmapCLParser::getGrowSize(CLParserPtr parser)
{
    return parser->get<size_t>("mmapGrowSize");
}
