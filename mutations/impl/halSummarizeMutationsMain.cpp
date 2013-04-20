/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halSummarizeMutations.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addOption("rootGenome", 
                           "name of root genome (none if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOption("maxGap", 
                           "maximum indel length to be considered a gap.  Gaps "
                           " can be nested within other rearrangements.", 
                           20);
  optionsParser->addOption("maxNFraction",
                           "maximum fraction of Ns in a rearranged segment "
                           "for it to not be ignored as missing data.",
                           1.0);
  optionsParser->addOptionFlag("justSubs", "just count substitutions. "
                               " Note results are total subs between genome "
                               " and all children, rather than branch results "
                               " when using the normal interface.  For tuning "
                               " and performance checking only", false);
  optionsParser->setDescription("Print summary table of mutation events "
                                "in the alignemt.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string rootGenomeName;
  string targetGenomes;
  string refSequenceName;
  hal_size_t maxGap;
  double nThreshold;
  bool justSubs;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    rootGenomeName = optionsParser->getOption<string>("rootGenome");
    targetGenomes = optionsParser->getOption<string>("targetGenomes");
    maxGap = optionsParser->getOption<hal_size_t>("maxGap");
    nThreshold = optionsParser->getOption<double>("maxNFraction");
    justSubs = optionsParser->getFlag("justSubs");

    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
                          " mutually exclusive");
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
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, 
                                                           optionsParser);
    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("hal alignmenet is empty");
    }
    
    set<const Genome*> targetSet;
    const Genome* rootGenome = NULL;
    if (rootGenomeName != "\"\"")
    {
      rootGenome = alignment->openGenome(rootGenomeName);
      if (rootGenome == NULL)
      {
        throw hal_exception(string("Root genome, ") + rootGenomeName + 
                            ", not found in alignment");
      }
      if (rootGenomeName != alignment->getRootName())
      {
        getGenomesInSubTree(rootGenome, targetSet);
      }
      else
      {
        rootGenome = NULL;
      }
    }

    if (targetGenomes != "\"\"")
    {
      vector<string> targetNames = chopString(targetGenomes, ",");
      for (size_t i = 0; i < targetNames.size(); ++i)
      {
        const Genome* tgtGenome = alignment->openGenome(targetNames[i]);
        if (tgtGenome == NULL)
        {
          throw hal_exception(string("Target genome, ") + targetNames[i] + 
                              ", not found in alignment");
        }
        targetSet.insert(tgtGenome);
      }
    }
    
    // convert back to strings since mtuations class will close genomes
    set<string> targetNames;
    if (targetSet.empty() == false)
    {
      for (set<const Genome*>::iterator i = targetSet.begin();
           i != targetSet.end(); ++i)
      {
        targetNames.insert((*i)->getName());
      }
    }
    
    SummarizeMutations mutations;
    mutations.analyzeAlignment(alignment, maxGap, nThreshold, justSubs,
                               targetSet.empty() ? NULL : &targetNames);

    cout << endl << mutations;
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
