/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "hal.h"

using namespace std;
using namespace hal;

static void extractAlignedRegions(const Genome* genome,
                                  ostream* bedStream);

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(false);
  optionsParser->setDescription("Extract aligned regions of genome (with its"
                                "parent) and output to bed file");
  optionsParser->addArgument("halPath", "input hal file");
  optionsParser->addArgument("genome", "name of genome to scan"
                             " (can't be root)");
  optionsParser->addOption("alignedFile", "path to bed file to write to", 
                           "stdout");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string genomeName;
  string outBedPath;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halPath");
    genomeName = optionsParser->getArgument<string>("genome");
    outBedPath = optionsParser->getOption<string>("alignedFile");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  try
  {
    AlignmentConstPtr inAlignment = openHalAlignmentReadOnly(halPath, 
                                                             optionsParser);
    if (inAlignment->getNumGenomes() == 0)
    {
      throw hal_exception("input hal alignmenet is empty");
    }
    
    const Genome* genome = inAlignment->openGenome(genomeName);
    
    if (genome == NULL)
    {
      throw hal_exception(string("Unable to open genome ") + genomeName +
                          "in alignment.");
    }
    if (genome->getParent() != NULL)
    {
      ostream* bedStream = &cout;
      bool newBed = false;
      if (outBedPath != "stdout")
      {
        bedStream = new ofstream(outBedPath.c_str());
        newBed = true;
      }
      if (!bedStream)
      {
        throw hal_exception(string("Error opening ") + outBedPath + 
                            " for writing");
      }
      extractAlignedRegions(genome, bedStream);
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

void extractAlignedRegions(const Genome* genome, ostream* bedStream)
{
  assert(genome && bedStream);
  if (genome->getNumTopSegments() > 0)
  {
    TopSegmentIteratorConstPtr topSeg = genome->getTopSegmentIterator();
    TopSegmentIteratorConstPtr endSeg = genome->getTopSegmentEndIterator();
    for (; topSeg->equals(endSeg) == false; topSeg->toRight())
    {
      if (topSeg->hasParent())
      {
        const Sequence* sequence = topSeg->getSequence();
        *bedStream << sequence->getName() << '\t'
                   << topSeg->getStartPosition() - sequence->getStartPosition()
                   << '\t' 
                   << topSeg->getEndPosition() - sequence->getStartPosition() +1
                   << '\n';
      }
    }
  }
  *bedStream << endl;
}
