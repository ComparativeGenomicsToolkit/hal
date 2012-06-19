/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halMafExport.h"

using namespace std;
using namespace hal;


int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("mafFile", "output maf file");
  optionsParser->addOption("refName", 
                           "name of reference genome (root if empty)", 
                           "\"\"");
  optionsParser->setDescription("Convert hal database to maf");
  string halPath;
  string mafPath;
  string refName;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    mafPath = optionsParser->getArgument<string>("mafFile");
    refName = optionsParser->getOption<string>("refName");
    if (refName == "\"\"")
    {
      refName = "";
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
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath);
    alignment->setOptionsFromParser(optionsParser);
    if (refName.empty())
    {
      refName = alignment->getRootName();
    }
    ofstream mafStream(mafPath.c_str());
    if (!mafStream)
    {
      throw hal_exception("Error opening " + mafPath);
    }
    MafExport mafExport;
    mafExport.convert(alignment, mafStream);
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
