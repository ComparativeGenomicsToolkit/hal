/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "hal4dExtract.h"


using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(false);
  optionsParser->setDescription("Extract Fourfold-Degenerate codon positions"
                                " from a BED file that contains exons");
  optionsParser->addArgument("halPath", "input hal file");
  optionsParser->addArgument("refGenome", "name of reference genome");
  optionsParser->addArgument("inBed", "path to bed file containing coding"
                             " exons in refGenome (or \"stdin\" to pipe from"
                             " standard input)");
  optionsParser->addArgument("outBed", "output path for bed file that "
                            "will only contain 4d sites (or \"stdout\" to "
                            "pipe to standard output)");
  optionsParser->addOption("bedVersion", "version of input bed file.  will be"
                           " automatically detected if not specified",
                           -1);
  optionsParser->addOptionFlag("append",
                               "append to instead of overwrite output file.",
                               false);
  optionsParser->addOptionFlag("conserved",
                               "ensure 4d sites are 4d sites in all leaf genomes",
                               false);
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string genomeName;
  string inBedPath;
  string outBedPath;
  int bedVersion;
  bool append;
  bool conserved;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halPath");
    genomeName = optionsParser->getArgument<string>("refGenome");
    inBedPath = optionsParser->getArgument<string>("inBed");
    outBedPath = optionsParser->getArgument<string>("outBed");
    bedVersion = optionsParser->getOption<int>("bedVersion");
    append = optionsParser->getFlag("append");
    conserved = optionsParser->getFlag("conserved");
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

    ifstream inBedFileStream;
    istream* inBedStream;
    if (inBedPath == "stdin")
    {
      inBedStream = &cin;
    }
    else
    {
      inBedFileStream.open(inBedPath.c_str());
      inBedStream = &inBedFileStream;
    }
    if (!(*inBedStream))
    {
      throw hal_exception("Error opening input " + inBedPath);
    }

    ofstream outBedFileStream;
    ostream* outBedStream;
    if (outBedPath == "stdout")
    {
      outBedStream = &cout;
    }
    else
    {
      ios_base::openmode openFlags = ios_base::out;
      if (append == true)
      {
        openFlags |= ios_base::app;
      }
      outBedFileStream.open(outBedPath.c_str(), openFlags);
      outBedStream = &outBedFileStream;
    }
    if (!(*outBedStream))
    {
      throw hal_exception("Error opening output " + inBedPath);
    }

    Extract4d extractor;
    extractor.run(genome, inBedStream, outBedStream, bedVersion, conserved);
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
