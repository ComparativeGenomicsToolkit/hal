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
                           "rearrangements in reference genome coordinates "
                           "(or stdout)",
                           "\"\"");
  optionsParser->addOption("parentFile", 
                           "bed file to write rearrangements "
                           "(deletions and duplications) in "
                           "reference's parent genome coordinates (or stdout)",
                           "\"\"");
  optionsParser->addOption("snpFile", 
                           "bed file write point mutations to "
                           "in reference genome coordinates (or stdout)",
                           "\"\"");
  optionsParser->addOption("delBreakFile", 
                           "bed file write deletion breakpoints to "
                           "in reference genome coordinates (or stdout)",
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
                           "maximum fraction of Ns in a rearranged segment "
                           "for it to not be ignored as missing data.",
                           1.0);
                           
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
    if (refBedPath == "\"\"" && parentBedPath == "\"\"" 
        && snpBedPath == "\"\"" && delBreakBedPath == "\"\"")
    {
      throw hal_exception("at least one of --refFile, --parentFile, "
                          "--snpFile or --delBreakFile must be specified");
    }
    if (refTargetsPath != "\"\"" && (
          refTargetsPath == refBedPath ||
          refTargetsPath == snpBedPath ||
          refTargetsPath == parentBedPath ||
          refTargetsPath == delBreakBedPath))
    {
      throw hal_exception("cannot output to same file as --refTargets");
    }
    if (parentBedPath != "\"\"" && (
          parentBedPath == refBedPath ||
          parentBedPath == snpBedPath ||
          parentBedPath == delBreakBedPath))
    {
      throw hal_exception("--parentBedPath must be unique");
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

    ostream* refBedStream = NULL;
    bool newRef = false;
    if (refBedPath != "\"\"")
    {
      if (refBedPath == "stdout")
      {
        refBedStream = &cout;
      }
      else
      {
        newRef = true;
        refBedStream = new ofstream(refBedPath.c_str());
      }
      if (!*refBedStream)
      {
        throw hal_exception("Error opening " + refBedPath);
      }
    }
  
    ostream* parentBedStream = NULL;
    bool newParent = false;
    if (parentBedPath != "\"\"")
    {
      if (parentBedPath == "stdout")
      {
        parentBedStream = &cout;
      }
      else
      {
        newParent = true;
        parentBedStream = new ofstream(parentBedPath.c_str());
      }
      if (!*parentBedStream)
      {
        throw hal_exception("Error opening " + parentBedPath);
      }
    }

    ostream* snpBedStream = NULL;
    bool newSnp = false;
    if (snpBedPath != "\"\"")
    {
      if (snpBedPath == refBedPath)
      {
        snpBedStream = refBedStream;
      }
      else if (snpBedPath == "stdout")
      {
        snpBedStream = &cout;
      }
      else
      {
        snpBedStream = new ofstream(snpBedPath.c_str());
        newSnp = true;
      }
      if (!*snpBedStream)
      {
        throw hal_exception("Error opening " + snpBedPath);
      }
    }

    ostream* delBreakBedStream = NULL;
    bool newDelBreak = false;
    if (delBreakBedPath != "\"\"")
    {
      if (delBreakBedPath == snpBedPath)
      {
        delBreakBedStream = snpBedStream;
      }
      else if (delBreakBedPath == refBedPath)
      {
        delBreakBedStream = refBedStream;
      }
      else if (delBreakBedPath == "stdout")
      {
        delBreakBedStream = &cout;
      }
      else
      {
        delBreakBedStream = new ofstream(delBreakBedPath.c_str());
        newDelBreak = true;
      }
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
        start = -1;
        end = (hal_size_t)-1;
        refSequenceName.erase();
        ss >> refSequenceName >> start >> end;
        if (ss.bad() || ss.fail() ||
            refSequenceName.empty() || start == -1 || end == (hal_size_t)-1)
        {
          if (!refSequenceName.empty() && refSequenceName[0] != '#')
          {
            cerr << "skipping malformed line " << line << " in " 
                 << refTargetsPath << "\n";
          }
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
                                    refBedStream, parentBedStream,
                                    snpBedStream, delBreakBedStream,
                                    refGenome, start, length);

          }
        }
      }
    }
    else
    {
      BranchMutations mutations;
      mutations.analyzeBranch(alignment, maxGap, nThreshold,
                              refBedStream, parentBedStream,
                              snpBedStream, delBreakBedStream,
                              refGenome, start, length);
    }
    if (newRef)
    {
      delete refBedStream;
    }
    if (newParent)
    {
      delete parentBedStream;
    }
    if (newSnp)
    {
      delete snpBedStream;
    }
    if (newDelBreak)
    {
      delete delBreakBedStream;
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
