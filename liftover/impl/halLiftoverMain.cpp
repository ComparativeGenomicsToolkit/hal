/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halLiftover.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("srcGenome", "source genome name");
  optionsParser->addArgument("srcBed", "path of input bed file");
  optionsParser->addArgument("tgtGenome", "target genome name");
  optionsParser->addArgument("tgtBed", "path of output bed file");
  optionsParser->addOptionFlag("addDupeColumn", "add column to output bed "
                               "with integer number of paralgous mappings in "
                               "target", false);
  optionsParser->addOptionFlag("append", "append results to tgtBed", false);
optionsParser->setDescription("Map genome interval coordinates between "
                              "different genomes.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string srcGenomeName;
  string srcBedPath;
  string tgtGenomeName;
  string tgtBedPath;
  bool dupeCol;
  bool append;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    srcGenomeName = optionsParser->getArgument<string>("srcGenome");
    srcBedPath =  optionsParser->getArgument<string>("srcBed");
    tgtGenomeName = optionsParser->getArgument<string>("tgtGenome");
    tgtBedPath =  optionsParser->getArgument<string>("tgtBed");
    dupeCol = optionsParser->getFlag("addDupeColumn");
    append = optionsParser->getFlag("append");
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
    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("hal alignmenet is empty");
    }

    const Genome* srcGenome = alignment->openGenome(srcGenomeName);
    if (srcGenome == NULL)
    {
      throw hal_exception(string("srcGenome, ") + srcGenomeName + 
                          ", not found in alignment");
    }
    const Genome* tgtGenome = alignment->openGenome(tgtGenomeName);
    if (tgtGenome == NULL)
    {
      throw hal_exception(string("tgtGenome, ") + tgtGenomeName + 
                          ", not found in alignment");
    }
    
    ifstream srcBed(srcBedPath.c_str());
    if (!srcBed)
    {
      throw hal_exception("Error opening srcBed, " + srcBedPath);
    }
    
    ios_base::openmode mode = append ? ios::out | ios::app : ios_base::out;
    ofstream tgtBed(tgtBedPath.c_str(), mode);
    if (!tgtBed)
    {
      throw hal_exception("Error opening tgtBed, " + tgtBedPath);
    }
    
    Liftover liftover;
    liftover.convert(alignment, srcGenome, &srcBed, tgtGenome, &tgtBed, 
                     dupeCol);

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
