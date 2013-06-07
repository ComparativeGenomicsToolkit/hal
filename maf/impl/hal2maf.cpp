/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "halMafExport.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("mafFile", "output maf file");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (root if empty)", 
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
                           "\"\"");
  optionsParser->addOption("refTargets", 
                           "bed file coordinates of intervals in the reference "
                           "genome to export",
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
  optionsParser->addOption("rootGenome", 
                           "name of root genome (none if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOption("maxRefGap", 
                           "maximum gap length in reference", 
                           0);
  optionsParser->addOptionFlag("noDupes", 
                               "ignore paralogy edges", 
                               false);
  optionsParser->addOptionFlag("noAncestors", 
                               "don't write ancestral sequences", 
                               false);
  optionsParser->addOptionFlag("ucscNames",
                               "use UCSC convention of Genome.Seqeunce "
                               "for output names.  By default, only sequence "
                               "names are used",
                               false);
  optionsParser->addOptionFlag("unique",
                               "only write column whose left-most reference "
                               "coordinate is in the specified range.  this "
                               "is used to insure that the same column isnt "
                               "sampled twice (due to ducplications) by mafs "
                               "generated on distinct ranges.",
                               false);
  optionsParser->addOptionFlag("append",
                               "append to instead of overwrite output file.",
                               false);

  optionsParser->setDescription("Convert hal database to maf.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string halPath;
  string mafPath;
  string refGenomeName;
  string rootGenomeName;
  string targetGenomes;
  string refSequenceName;
  string refTargetsPath;
  hal_index_t start;
  hal_size_t length;
  hal_size_t maxRefGap;
  bool noDupes;
  bool noAncestors;
  bool ucscNames;
  bool unique;
  bool append;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    mafPath = optionsParser->getArgument<string>("mafFile");
    refGenomeName = optionsParser->getOption<string>("refGenome");
    rootGenomeName = optionsParser->getOption<string>("rootGenome");
    targetGenomes = optionsParser->getOption<string>("targetGenomes");
    refSequenceName = optionsParser->getOption<string>("refSequence");    
    refTargetsPath = optionsParser->getOption<string>("refTargets");
    start = optionsParser->getOption<hal_index_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    maxRefGap = optionsParser->getOption<hal_size_t>("maxRefGap");
    noDupes = optionsParser->getFlag("noDupes");
    noAncestors = optionsParser->getFlag("noAncestors");
    ucscNames = optionsParser->getFlag("ucscNames");
    unique = optionsParser->getFlag("unique");
    append = optionsParser->getFlag("append");

    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
                          "mutually exclusive");
    }
    if (refTargetsPath != "\"\"" && refSequenceName != "\"\"")
    {
      throw hal_exception("--refSequence and --refTargets options are "
                          "mutually exclusive");
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

    const Genome* refGenome = NULL;
    if (refGenomeName != "\"\"")
    {
      refGenome = alignment->openGenome(refGenomeName);
      if (refGenome == NULL)
      {
        throw hal_exception(string("Reference genome, ") + refGenomeName + 
                            ", not found in alignment");
      }
    }
    else
    {
      refGenome = alignment->openGenome(alignment->getRootName());
    }
    const SegmentedSequence* ref = refGenome;
    
    const Sequence* refSequence = NULL;
    if (refSequenceName != "\"\"")
    {
      refSequence = refGenome->getSequence(refSequenceName);
      ref = refSequence;
      if (refSequence == NULL)
      {
        throw hal_exception(string("Reference sequence, ") + refSequenceName + 
                            ", not found in reference genome, " + 
                            refGenome->getName());
      }
    }

    ios_base::openmode openFlags = ios_base::out;
    if (append == true)
    {
      openFlags |= ios_base::app;
    }
    ofstream mafStream(mafPath.c_str(), openFlags);
    if (!mafStream)
    {
      throw hal_exception("Error opening " + mafPath);
    }

    MafExport mafExport;
    mafExport.setMaxRefGap(maxRefGap);
    mafExport.setNoDupes(noDupes);
    mafExport.setNoAncestors(noAncestors);
    mafExport.setUcscNames(ucscNames);
    mafExport.setUnique(unique);

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
            refSequenceName.empty() || start == -1 || end == (hal_size_t)-1 ||
            (hal_index_t)end <= start)
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
          if (refSequence != NULL && 
              start + length <= refSequence->getSequenceLength())
          {
            mafExport.convertSegmentedSequence(mafStream, alignment, 
                                               refSequence, 
                                               start, length, targetSet);
          }
          else
          {
            cerr << "bed coordinate " << refSequenceName << " " 
                 << start << " " << end << " not found in " 
                 << refGenome->getName() << "\n";
          }
        }
      }
    }
    else
    {
      mafExport.convertSegmentedSequence(mafStream, alignment, ref, 
                                         start, length, targetSet);
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
