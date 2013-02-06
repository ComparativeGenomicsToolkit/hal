/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include "halChain.h"

using namespace std;
using namespace hal;

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->setDescription("Rertrieve chain (pairwise alignment) "
                                "information from a hal database");
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addArgument("genome", "(query) genome to process");
  optionsParser->addOption("sequence", "sequence name in query genome ("
                           "all sequences if not specified)", "\"\"");
  optionsParser->addOption("start", "start position in query genome", 0);
  optionsParser->addOption("length", "maximum length of chain to output.", 0);
  optionsParser->addOption("chainFile", "path for output file.  stdout if not"
                           " specified", "\"\"");

  string halPath;
  string chainPath;
  string genomeName;
  string sequenceName;
  hal_size_t start;
  hal_size_t length;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    genomeName = optionsParser->getArgument<string>("genome");
    sequenceName = optionsParser->getOption<string>("sequence");
    start = optionsParser->getOption<hal_size_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    chainPath = optionsParser->getOption<string>("chainFile");
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
      throw hal_exception(string("Genome not found: ") + genomeName);
    }
    
    const Sequence* sequence = NULL;
    if (sequenceName != "\"\"")
    {
      sequence = genome->getSequence(sequenceName);
      if (sequence == NULL)
      {
        throw hal_exception(string("Sequence not found: ") + sequenceName);
      }
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


