/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "hal.h"
#include "halChain.h"
#include "halBlockInterpolate.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halPath", "input hal file");
  optionsParser->addArgument("genome", "reference genome to scan");
  optionsParser->addOption("sequence", "sequence name to scan ("
                           "all sequences by default)", 
                           "\"\"");
  optionsParser->addOption("start",
                           "coordinate within reference genome (or sequence"
                           " if specified) to start at",
                           0);
  optionsParser->addOption("length",
                           "length of the reference genome (or sequence"
                           " if specified) to convert.  If set to 0,"
                           " the entire thing is converted",
                           0);
  optionsParser->addOption("step", "step size", 1);
  optionsParser->addOption("targetGenome", "target genome [default: "
                           "parent of reference genome", "\"\"");
  
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string genomeName;
  string tgtGenomeName;
  string sequenceName;
  hal_size_t start;
  hal_size_t length;
  hal_size_t step;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halPath");
    genomeName = optionsParser->getArgument<string>("genome");
    sequenceName = optionsParser->getOption<string>("sequence");
    start = optionsParser->getOption<hal_size_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    step = optionsParser->getOption<hal_size_t>("step");
    tgtGenomeName = optionsParser->getOption<string>("targetGenome");
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
      throw hal_exception("input hal alignmenet is empty");
    }
    
    const Genome* genome = alignment->openGenome(genomeName);
    if (genome == NULL)
    {
      throw hal_exception(string("Genome ") + genomeName + " not found");
    }

    const Genome* target;
    if (tgtGenomeName != "\"\"")
    {
      target = alignment->openGenome(tgtGenomeName);
      if (target == NULL)
      {
        throw hal_exception(string("Genome ") + tgtGenomeName + " not found");
      }
    }
    else
    {
      target = genome->getParent();
      if (target == NULL)
      {
        throw hal_exception(string("Genome ") + genomeName + " has no parent: "
                            "--targetGenome must be specified");
      }
    }

    const Sequence* sequence = NULL;
    if (sequenceName != "\"\"")
    {
      sequence = genome->getSequence(sequenceName);
      if (sequence == NULL)
      {
        throw hal_exception(string("Sequence ") + sequenceName + " not found");
      }
    }
    
    blockInterpolateGenome(genome, sequence, target, start, length, step);
    
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

