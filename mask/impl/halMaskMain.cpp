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
  optionsParser->setDescription("Rertrieve basic statics from a hal database");
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addArgument("maskBedFile", "path to bed file to write to");
  optionsParser->addOption("extend", "extend masked regions by given num. "
                           "of bases.", 0);
  optionsParser->addOption("extendPct", "extend masked regions by percentage"
                           " of their lengths", 0);

  string halPath;
  string bedPath;
  hal_size_t extend;
  double extendPct;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    bedPath = optionsParser->getArgument<string>("maskBedFile");
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

    ofstream bedStream(bedPath.c_str());
    if (!bedStream)
    {
      throw hal_exception(string("Error opening ") + bedPath + " for writing");
    }

    MaskExtractor mask;
    mask.extract(alignment, &bedStream, extend, extendPct);
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

