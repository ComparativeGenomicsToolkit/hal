/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halMafScanDimensions.h"
#include "halMafWriteGenomes.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("mafFile", "output maf file");
  optionsParser->addArgument("refGenome", "name of reference genome in MAF");
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOptionFlag("append",
                               "append maf as subtree to existing alignment."
                               " reference must alaready be present in hal"
                               " dabase as a leaf.",
                               false);
                           
  optionsParser->setDescription("import maf into hal database.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string halPath;
  string mafPath;
  string refGenomeName;
  string targetGenomes;
  bool append;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    mafPath = optionsParser->getArgument<string>("mafFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    targetGenomes = optionsParser->getOption<string>("targetGenomes");
    append = optionsParser->getFlag("append");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }
//  try
  {
    ifstream mafStream(mafPath.c_str());
    if (!mafStream)
    {
      throw hal_exception("Error opening MAF file: " + mafPath);
    }    
    if (!ofstream(halPath.c_str(), ios_base::app))
    {
      throw hal_exception("Error opening HAL file: " + halPath);
    }
    
    AlignmentPtr alignment;
    if (append == true)
    {
      alignment = openHalAlignment(halPath, optionsParser);
      if (alignment->openGenome(refGenomeName) == NULL)
      {
        throw hal_exception("Reference genome " + refGenomeName + " not found "
                            "in hal file");
      }
      else if (alignment->openGenome(refGenomeName)->getNumChildren() > 0)
      {
        throw hal_exception("Reference genome " + refGenomeName + " not a leaf "
                            "in hal file");
      }
    }
    else
    {
      alignment = hdf5AlignmentInstance();
      alignment->setOptionsFromParser(optionsParser);
      alignment->createNew(halPath);
    }
    
    vector<string> targetNames;
    if (targetGenomes != "\"\"")
    {
      targetNames = chopString(targetGenomes, ",");
    }
    set<string> targetSet(targetNames.begin(), targetNames.end());
    targetSet.insert(refGenomeName);

    MafScanDimensions dScan;
    dScan.scan(mafPath, targetSet);
    MafWriteGenomes writer;
    writer.convert(mafPath, refGenomeName, targetSet, dScan.getDimensions(),
                   alignment);
    
    const MafScanDimensions::DimMap& dimMap = dScan.getDimensions();
    for (MafScanDimensions::DimMap::const_iterator i = dimMap.begin();
         i != dimMap.end(); ++i)
    {
      cout << i->first << ":  " << i->second->_length << ", "
           << i->second->_segments << ", " << i->second->_intervals.size()
           << "\n";
    }

  }try{}
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
