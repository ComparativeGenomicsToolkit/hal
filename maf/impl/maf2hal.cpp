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
#include "halMafScanReference.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("mafFile", "output maf file");
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addOption("refGenome", "name of reference genome in MAF "
                           "(first found if empty)",
                           "\"\"");
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
    refGenomeName = optionsParser->getOption<string>("refGenome");
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

    if (refGenomeName ==  "\"\"")
    {
      MafScanReference nameGetter;
      refGenomeName = nameGetter.getRefName(mafPath);
      if (refGenomeName.empty())
      {
        throw hal_exception("Unable to find a single genome name in maf file");
      }
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

    string prevGenome, curGenome;
    hal_size_t segmentCount = 0;
    hal_size_t setCount = 0;
    hal_size_t sequenceCount = 0;
    hal_size_t skipCount = 0;
    const MafScanDimensions::DimMap& dimMap = dScan.getDimensions();
    for (MafScanDimensions::DimMap::const_iterator i = dimMap.begin();
         i != dimMap.end(); ++i)
    {
      curGenome = MafScanner::genomeName(i->first);
      MafScanDimensions::DimMap::const_iterator next = i;
      ++next;
      if (prevGenome.empty() == false && (curGenome != prevGenome ||
                                          next == dimMap.end()))
      {
        cout << prevGenome << ":  " <<  segmentCount << " segments, "
             << setCount << " set entries, " << skipCount << " dupe rows, "
             << sequenceCount << " sequences"
             << endl;
        segmentCount = 0;
        setCount = 0;
        skipCount = 0;
        sequenceCount = 0;
      }
      segmentCount += i->second->_numSegments;
      setCount += i->second->_startMap.size();
      sequenceCount++;
      skipCount += i->second->_badPosSet.size();
      prevGenome = curGenome;
    }
    cout << "Total Number of blocks in maf: " << dScan.getNumBlocks() << "\n";

    MafWriteGenomes writer;
    writer.convert(mafPath, refGenomeName, targetSet, dScan.getDimensions(),
                   alignment);


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
