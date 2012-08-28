/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include "halCons.h"

using namespace std;
using namespace hal;


int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addOption("gapThreshold", "max len of gap indel",
    10);
  optionsParser->setDescription("Gather preliminary conservation statistics"
                                " from hal database");
  string path;
  hal_size_t gapThreshold;
  try
  {
    optionsParser->parseOptions(argc, argv);
    path = optionsParser->getArgument<string>("halFile");
    gapThreshold = optionsParser->getOption<hal_size_t>("gapThreshold");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }
  try
  {
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(path);
    alignment->setOptionsFromParser(optionsParser);
    HalCons halCons(alignment, gapThreshold);
    cout << endl << halCons;
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
