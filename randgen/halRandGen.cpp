/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctime>

#include <H5Cpp.h>
#include "hal.h"
#include "halRandomData.h"


using namespace std;
using namespace hal;
using namespace H5;

struct Options {
   double _meanDegree;
   double _maxBranchLength;
   hal_size_t _maxGenomes;
   hal_size_t _minSegmentLength;
   hal_size_t _maxSegmentLength;
   hal_size_t _minSegments;
   hal_size_t _maxSegments;
   hsize_t _hdf5Chunk;
   hsize_t _hdf5Compression;
   int _seed;
};

static const Options defaultSm = { 0.75, 0.1, 5, 10, 1000, 5, 10, 2000000, 9, -1};
static const Options defaultMed = { 1.25,  0.7, 20, 2, 50, 1000, 50000, 
                                    2000000, 9, -1 };
static const Options defaultBig = { 2,  0.7, 50, 2, 500, 100, 5000, 
                                    2000000, 9, -1 };
static const Options defaultLrg = {2, 1, 100, 2, 10, 10000, 500000,
                                   2000000, 9, -1};

static void printUsage()
{
  cerr << "Usage: halStats [options] <path of ouput hal alignment file>\n"
       << "[options]:\n"
       << "--preset <small, medium, big, large> [medum]\n"
       << "--meanDegree <double> [" << defaultMed._meanDegree << "]\n"
       << "--maxBranchLength <double> [" << defaultMed._maxBranchLength 
       << "]\n"
       << "--maxGenomes <int> [" << defaultMed._maxGenomes << "]\n"
       << "--minSegmentLength <int> [" << defaultMed._minSegmentLength << "]\n"
       << "--maxSegmentLength <int> [" << defaultMed._maxSegmentLength << "]\n"
       << "--maxSegments <int> [" << defaultMed._maxSegments << "]\n"
       << "--minSegments <int> [" << defaultMed._minSegments << "]\n"
       << "--hdf5Chunk <int> [" << defaultMed._hdf5Chunk << "]\n"
       << "--hdf5Compression <int> [" << defaultMed._hdf5Compression << "]\n"
       << "--seed <int> [system time]\n"
       << endl;
  exit(1);
}

static Options checkOptions(int argc, char** argv)
{
  if (argc == 1 || string(argv[1]) == "--help")
  {
    printUsage();
  }

  Options options = defaultMed;

  for (int i = 1; i < argc - 1; ++i)
  {
    if (i == argc - 2) 
    {
      printUsage();
    }

    string arg(argv[i]);
    stringstream val;
    val << argv[++i];
    
    if (arg == "--preset")
    {
      string preset;
      val >> preset;
      if (preset == "small") options = defaultSm;
      else if (preset == "medium") options = defaultMed;
      else if (preset == "big") options = defaultBig;
      else if (preset == "large") options = defaultLrg;
      else printUsage();
    }
    else if (arg == "--meanDegree")
    {
      val >> options._meanDegree; 
    }
    else if (arg == "--maxBranchLength")
    {
      val >> options._maxBranchLength;
    }
    else if (arg == "--maxGenoems")
    {
      val >> options._maxGenomes;
    }
    else if (arg == "--minSegmentLength")
    {
      val >> options._minSegmentLength;
    }
    else if (arg == "--maxSegmentLength")
    {
      val >> options._maxSegmentLength;
    }
    else if (arg == "--minSegments")
    {
      val >> options._minSegments;
    }
    else if (arg == "--maxSegments")
    {
      val >> options._maxSegments;
    }
    else if (arg == "--hdf5Chunk")
    {
      val >> options._hdf5Chunk;
    }
    else if (arg == "--hdf5Compression")
    {
      val >> options._hdf5Compression;
    }
    else if (arg == "--seed")
    {
      val >> options._seed;
    }

    else
    {
      printUsage();
    }
  }
  return options;
}
 
int main(int argc, char** argv)
{
  Options options;
  string outPath;
  try
  {
    options = checkOptions(argc, argv);
    outPath = argv[argc - 1];
  }
  catch (...)
  {
    printUsage();
  }

  try
  {
    // load up the hdf5 options
    DSetCreatPropList dcprops = DSetCreatPropList::DEFAULT;
    if (options._hdf5Compression < 10)
    {
      dcprops.setDeflate(options._hdf5Compression);      
    }
    // can't get szip to create.  must be violating some rule but i don't know
    // what.  should explore when i have time.  
    else if (options._hdf5Compression < 200)
    {
      dcprops.setSzip(H5_SZIP_NN_OPTION_MASK,32);
    }
    else
    {
      dcprops.setSzip(H5_SZIP_NN_OPTION_MASK, options._hdf5Compression - 200);
    }

    dcprops.setChunk(1, &options._hdf5Chunk);

    FileAccPropList aprops = FileAccPropList::DEFAULT;
    aprops.setCache(11, 5, 10000000, 0.25);

    AlignmentPtr alignment = hdf5AlignmentInstance(FileCreatPropList::DEFAULT,
                                                   aprops,
                                                   dcprops);

    alignment->createNew(outPath);

    // call the crappy unit-test simulator 
    createRandomAlignment(alignment,
                          options._meanDegree,
                          options._maxBranchLength,
                          options._maxGenomes,
                          options._minSegmentLength,
                          options._maxSegmentLength,
                          options._minSegments,
                          options._maxSegments,
                          options._seed);
    
    alignment->close();
  }
  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
  
  return 0;
}
