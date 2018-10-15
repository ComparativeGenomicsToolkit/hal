/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "mmapCLParser.h"
#include "mmapFile.h"

using namespace hal;
using namespace std;

void MMapCLParser::defineOptions(CLParserPtr parser,
                                 unsigned mode)
{
    if (mode & CREATE_ACCESS) {
      parser->addOption("mmapInitSize", "mmap HAL file initial size", MMAP_DEFAULT_INIT_SIZE);
    }
    if (mode & (CREATE_ACCESS | WRITE_ACCESS)) {
      parser->addOption("mmapGrowSize", "mmap HAL file size to grow when more memory is needed", MMAP_DEFAULT_INIT_SIZE);
    }
}

size_t MMapCLParser::getInitSize(CLParserPtr parser)
{
    if (parser->hasOption("mmapInitSize")) {
        return parser->get<size_t>("mmapInitSize");
    } else {
        return 0;
    }
}

size_t MMapCLParser::getGrowSize(CLParserPtr parser)
{
    if (parser->hasOption("mmapGrowSize")) {
        return parser->get<size_t>("mmapGrowSize");
    } else {
        return 0;
    }
}
