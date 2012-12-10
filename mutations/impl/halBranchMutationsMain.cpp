/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halBranchMutations.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("refGenome", 
                             "name of reference genome (analyzed branch is "
                             "this genome and its parent).");
  optionsParser->addOption("refFile", 
                           "bed file to write structural "
                           "rearrangements in reference genome coordinates",
                           "\"\"");
  optionsParser->addOption("parentFile", 
                           "bed file to write rearrangements "
                           "(deletions and duplications) in "
                           "reference's parent genome coordinates",
                           "\"\"");
  optionsParser->addOption("snpFile", 
                           "bed file write point mutations to "
                           "in reference genome coordinates",
                           "\"\"");
  optionsParser->addOption("delBreakFile", 
                           "bed file write deletion breakpoints to "
                           "in reference genome coordinates",
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
                           "\"\"");
  optionsParser->addOption("refTargets",
                           "bed file coordinates of intervals in the reference "
                           "genome to analyze",
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
  optionsParser->addOption("maxGap", 
                           "maximum indel length to be considered a gap.  Gaps "
                           " can be nested within other rearrangements.", 
                           20);
  optionsParser->addOption("maxNFraction",
                           "maximum franction of Ns in a rearranged segment "
                           "for it to not be ignored as missing data.",
                           .10);
                           
  optionsParser->setDescription("Identify mutations on branch between given "
                                "genome and its parent.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string refBedPath;
  string parentBedPath;
  string snpBedPath;
  string delBreakBedPath;
  string refGenomeName;
  string refSequenceName;
  string refTargetsPath;
  hal_index_t start;
  hal_size_t length;
  hal_size_t maxGap;
  double nThreshold;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    refBedPath = optionsParser->getOption<string>("refFile");
    parentBedPath = optionsParser->getOption<string>("parentFile");
    snpBedPath = optionsParser->getOption<string>("snpFile");
    delBreakBedPath = optionsParser->getOption<string>("delBreakFile");
    refSequenceName = optionsParser->getOption<string>("refSequence");
    refTargetsPath = optionsParser->getOption<string>("refTargets");
    start = optionsParser->getOption<hal_index_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    maxGap = optionsParser->getOption<hal_size_t>("maxGap");
    nThreshold = optionsParser->getOption<double>("maxNFraction");
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
    
    const Genome* refGenome = NULL;
    refGenome = alignment->openGenome(refGenomeName);
    if (refGenome == NULL)
    {
      throw hal_exception(string("Reference genome, ") + refGenomeName + 
                          ", not found in alignment");
    }
    if (refGenome->getName() == alignment->getRootName())
    {
      throw hal_exception("Reference genome must denote bottom node in "
                          "a branch, and therefore cannot be the root.");
    }
    if (start + length >= refGenome->getSequenceLength())
    {
      throw hal_exception(string("Invalid range for ") + refGenomeName);
    }
    if (refBedPath == "\"\"" && parentBedPath == "\"\"" && snpBedPath == "\"\"")
    {
      throw hal_exception("at least one of --refFile, --parentFile or "
                          "--snpFile must be specified");
    }

    const Sequence* refSequence = NULL;
    if (refSequenceName != "\"\"")
    {
      refSequence = refGenome->getSequence(refSequenceName);
      if (refSequence == NULL)
      {
        throw hal_exception(string("Reference sequence, ") + refSequenceName + 
                            ", not found in reference genome, " + 
                            refGenome->getName());
      }
      if (start + length >= refSequence->getSequenceLength())
      {
        throw hal_exception(string("Invalid range for ") + refSequenceName);
      }
      if (length == 0)
      {
        length = refSequence->getSequenceLength() - start;
      }
      start += refSequence->getStartPosition();
    }

    if (length == 0)
    {
      assert(refSequence == NULL);
      length = refGenome->getSequenceLength() - start;
    }

    ofstream refBedStream;
    if (refBedPath != "\"\"")
    {
      refBedStream.open(refBedPath.c_str());
      if (!refBedStream)
      {
        throw hal_exception("Error opening " + refBedPath);
      }
    }
  
    ofstream parentBedStream;
    if (parentBedPath != "\"\"")
    {
      parentBedStream.open(parentBedPath.c_str());
      if (!parentBedStream)
      {
        throw hal_exception("Error opening " + parentBedPath);
      }
    }

    ofstream snpBedStream;
    if (snpBedPath != "\"\"")
    {
      snpBedStream.open(snpBedPath.c_str());
      if (!snpBedStream)
      {
        throw hal_exception("Error opening " + snpBedPath);
      }
    }

    ofstream delBreakBedStream;
    if (delBreakBedPath != "\"\"")
    {
      delBreakBedStream.open(delBreakBedPath.c_str());
      if (!delBreakBedStream)
      {
        throw hal_exception("Error opening " + delBreakBedPath);
      }
    }

    ifstream refTargetsStream;
    if (refTargetsPath != "\"\"")
    {
      refTargetsStream.open(refTargetsPath.c_str());
      if (!refTargetsStream)
      {
        throw hal_exception("Error opening " + refTargetsPath);
      }
      string line;
      hal_size_t end;
      while (!refTargetsStream.bad() && !refTargetsStream.eof())
      {
        getline(refTargetsStream, line);
        istringstream ss(line);
        ss >> refSequenceName >> start >> end;
        if (ss.bad() && (!refSequenceName.empty() || refSequenceName[0] != '#'))
        {
          cerr << "skipping malformed line " << line << " in " 
               << refTargetsPath << "\n";
        }
        else
        {
          refSequence = refGenome->getSequence(refSequenceName);
          length = end - start;
          if (refSequence != NULL && length <= refSequence->getSequenceLength())
          {
            start = refSequence->getStartPosition() + start;
            BranchMutations mutations;
            mutations.analyzeBranch(alignment, maxGap, nThreshold,
                                    refBedStream.is_open() ? &refBedStream : 
                                    NULL, 
                                    parentBedStream.is_open() ? 
                                    &parentBedStream : NULL ,
                                    snpBedStream.is_open() ? &snpBedStream : 
                                    NULL, 
                                    delBreakBedStream.is_open() ? 
                                    &delBreakBedStream : NULL,
                                    refGenome, start, length);

          }
        }
      }
    }
    else
    {
      BranchMutations mutations;
      mutations.analyzeBranch(alignment, maxGap, nThreshold,
                              refBedStream.is_open() ? &refBedStream : NULL, 
                              parentBedStream.is_open() ? &parentBedStream : 
                              NULL ,
                              snpBedStream.is_open() ? &snpBedStream : NULL, 
                              delBreakBedStream.is_open() ? 
                              &delBreakBedStream : NULL,
                              refGenome, start, length);
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
