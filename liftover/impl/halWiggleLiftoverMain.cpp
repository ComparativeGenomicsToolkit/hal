/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halWiggleLiftover.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("srcGenome", "source genome name");
  optionsParser->addArgument("srcWig", "path of input .wig file.  set as stdin "
                             "to stream from standard input");
  optionsParser->addArgument("tgtGenome", "target genome name");
  optionsParser->addArgument("tgtWig", "path of output .wig file.  set as stdout"
                             " to stream to standard output.");
  optionsParser->addOptionFlag("noDupes", "do not map between duplications in"
                               " graph.", false);
  optionsParser->addOptionFlag("append", "append/merge results into tgtWig.  "
                               "Note that the entire tgtWig file will be loaded into"
                               " memory then overwritten, so this data can be lost "
                               "in event of a crash", false);
/*  optionsParser->addOptionFlag("unique",
                               "only map block if its left-most paralog is in"
                               "the input.  this "
                               "is used to insure that the same column isnt "
                               "sampled twice (due to ducplications) by mafs "
                               "generated on distinct ranges.",
                               false);
*/
  optionsParser->setDescription("Map wiggle genome annotation between two"
                                " genomes.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string srcGenomeName;
  string srcWigPath;
  string tgtGenomeName;
  string tgtWigPath;
  bool noDupes;
  bool append;
  bool unique;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    srcGenomeName = optionsParser->getArgument<string>("srcGenome");
    srcWigPath =  optionsParser->getArgument<string>("srcWig");
    tgtGenomeName = optionsParser->getArgument<string>("tgtGenome");
    tgtWigPath =  optionsParser->getArgument<string>("tgtWig");
    noDupes = optionsParser->getFlag("noDupes");
    append = optionsParser->getFlag("append");
    //  unique = optionsParser->getFlag("unique");
    unique = false;
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
      throw hal_exception("hal alignmnet is empty");
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
    
    ifstream srcWig;
    istream* srcWigPtr;
    if (srcWigPath == "stdin")
    {
      srcWigPtr = &cin;
    }
    else
    {
      srcWig.open(srcWigPath.c_str());
      srcWigPtr = &srcWig;
      if (!srcWig)
      {
        throw hal_exception("Error opening srcWig, " + srcWigPath);
      }
    }

    WiggleLiftover liftover;
    if (append == true && tgtWigPath != "stdout")
    {
      // load the wig data into memory so that it can be properly merged
      // with the new data from the liftover.
      ifstream tgtWig(tgtWigPath.c_str());
      if (tgtWig)
      {
        liftover.preloadOutput(alignment, tgtGenome, &tgtWig);
      }
    }

    ofstream tgtWig;
    ostream* tgtWigPtr;
    if (tgtWigPath == "stdout")
    {
      tgtWigPtr = &cout;
    }
    else
    {      
      tgtWig.open(tgtWigPath.c_str());
      tgtWigPtr = &tgtWig;
      if (!tgtWig)
      {
        throw hal_exception("Error opening tgtWig, " + tgtWigPath);
      }
    }

    liftover.convert(alignment, srcGenome, srcWigPtr, tgtGenome, tgtWigPtr,
                     !noDupes, unique);
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
