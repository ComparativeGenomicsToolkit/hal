/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halMaskExtractor.h"

using namespace std;
using namespace hal;

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->setDescription("Write masked intervals of genome into bed "
                                "file");
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addArgument("genome", "name of genome to process");
  optionsParser->addOption("maskFile", "path to bed file to write to", 
                           "stdout");
  optionsParser->addOption("extend", "extend masked regions by given num. "
                           "of bases.", 0);
  optionsParser->addOption("extendPct", "extend masked regions by percentage"
                           " of their lengths", 0);

  string halPath;
  string genomeName;
  string bedPath;
  hal_size_t extend;
  double extendPct;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    genomeName = optionsParser->getArgument<string>("genome");
    bedPath = optionsParser->getOption<string>("maskFile");
    extend = optionsParser->getOption<hal_size_t>("extend");
    extendPct = optionsParser->getOption<double>("extendPct");

    if (extend != 0 && extendPct != 0.)
    {
      throw hal_exception("--extend and --extendPct options are exclusive.");
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }
  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, 
                                                           optionsParser);

    const Genome* genome = alignment->openGenome(genomeName);
    if (genome == NULL)
    {
      throw hal_exception(string("Genome ") + genomeName + " not found.");
    }

    ostream* bedStream = &cout;
    bool newBed = false;
    if (bedPath != "stdout")
    {
      bedStream = new ofstream(bedPath.c_str());
      newBed = true;
    }
    if (!bedStream)
    {
      throw hal_exception(string("Error opening ") + bedPath + " for writing");
    }

    MaskExtractor mask;
    mask.extract(alignment, genome, bedStream, extend, extendPct);

    if (newBed)
    {
      delete bedStream;
    }
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

