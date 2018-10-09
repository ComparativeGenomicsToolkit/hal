/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include <cstdlib>
#include <deque>
#include "hdf5CLParser.h"

using namespace hal;
using namespace std;
using namespace H5;

const hsize_t HDF5CLParser::DefaultChunkSize = 1000;
const hsize_t HDF5CLParser::DefaultDeflate = 2;
const hsize_t HDF5CLParser::DefaultCacheMDCElems = 113;
const hsize_t HDF5CLParser::DefaultCacheRDCElems = 599999;
const hsize_t HDF5CLParser::DefaultCacheRDCBytes = 15728640;
const double HDF5CLParser::DefaultCacheW0 = 0.75;
const bool HDF5CLParser::DefaultInMemory = false;

void HDF5CLParser::defineOptions(CLParserPtr parser,
                                 bool createOptions)
{
  if (createOptions)
  {
    parser->addOption("hdf5Chunk", "hdf5 chunk size", DefaultChunkSize);
    parser->addOption("chunk", " obsolete name for --hdf5Chunk ", DefaultChunkSize);

    parser->addOption("hdf5Compression", "hdf5 compression factor [0:none - 9:max]", 
                      DefaultDeflate);
    parser->addOption("deflate", "obsolete name for --hdf5Compression", 
                      DefaultDeflate);
  }
  parser->addOption("hdf5CacheMDC", "number of metadata slots in hdf5 cache",
                      DefaultCacheMDCElems);
  parser->addOption("cacheMDC", "obsolete name for --hdf5CacheMDC ",
                      DefaultCacheMDCElems);

  parser->addOption("hdf5CacheRDC", "number of regular slots in hdf5 cache.  should be"
                      " a prime number ~= 10 * DefaultCacheRDCBytes / chunk",
                      DefaultCacheRDCElems);
  parser->addOption("cacheRDC", "obsolete name for --hdf5CacheRDC",
                      DefaultCacheRDCElems);

  parser->addOption("hdf5CacheBytes", "maximum size in bytes of regular hdf5 cache",
                      DefaultCacheRDCBytes);
  parser->addOption("cacheBytes", "obsolete name for --hdf5CacheBytes",
                      DefaultCacheRDCBytes);

  parser->addOption("hdf5CacheW0", "w0 parameter for hdf5 cache", DefaultCacheW0);
  parser->addOption("cacheW0", "obsolete name for --hdf5CacheW0", DefaultCacheW0);

  parser->addOptionFlag("hdf5InMemory", "load all data in memory (and disable hdf5 cache)",
                        DefaultInMemory);
  parser->addOptionFlag("inMemory", "obsolete name for --hdf5InMemory",
                        DefaultInMemory);
#ifdef ENABLE_UDC
  // this can be define by other storage formats as well
  if (not parser->hasArgument("udcCacheDir")) {
      parser->addOption("udcCacheDir", "udc cache path for *input* hal file(s).",
                        "\"\"");
  }
#endif
}

template <typename T>
hsize_t getOptionAlt(CLParserPtr parser,
                     const std::string& name,
                     const std::string& obsoleteName) {
    if (parser->specifiedOption(obsoleteName)) {
        cerr << "Warning: --" << obsoleteName << " is obsolete, use --" << name << endl;
        return parser->getOption<T>(obsoleteName);
    } else {
        return parser->getOption<T>(name);
    }
}

static bool getFlagAlt(CLParserPtr parser,
                       const std::string& name,
                       const std::string& obsoleteName) {
    if (parser->specifiedFlag(obsoleteName)) {
        cerr << "Warning: --" << obsoleteName << " is obsolete, use --" << name << endl;
        return parser->getFlag(obsoleteName);
    } else {
        return parser->getFlag(name);
    }
}


void HDF5CLParser::applyToDCProps(CLParserPtr parser,
                                  DSetCreatPropList& dcprops)
{
    hsize_t chunk = getOptionAlt<hsize_t>(parser, "hdf5Chunk", "chunk");
    hsize_t deflate = getOptionAlt<hsize_t>(parser, "hdf5Compression", "deflate");
    dcprops.setChunk(1, &chunk);
    dcprops.setDeflate(deflate);
}

void HDF5CLParser::applyToAProps(CLParserPtr parser,
                                 H5::FileAccPropList& aprops)
{
    hsize_t mdc = getOptionAlt<hsize_t>(parser, "hdf5CacheMDC", "cacheMDC");
    hsize_t rdc = getOptionAlt<hsize_t>(parser, "hdf5CacheRDC", "cacheRDC");
    hsize_t rdcb = getOptionAlt<hsize_t>(parser, "hdf5CacheBytes", "cacheBytes");
    double w0 = getOptionAlt<double>(parser, "hdf5CacheW0", "cacheW0");
    aprops.setCache(mdc, rdc, rdcb, w0);
}

bool HDF5CLParser::getInMemory(CLParserPtr parser)
{
    return getFlagAlt(parser, "hdf5InMemory", "inMemory");}
