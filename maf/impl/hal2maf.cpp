/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "halMafExport.h"
#include "halMafBed.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("mafFile", "output maf file (or \"stdout\" to "
                             "pipe to standard output)");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (root if empty)", 
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
                           "\"\"");
  optionsParser->addOption("refTargets", 
                           "bed file coordinates of intervals in the reference "
                           "genome to export (or \"stdin\" to pipe from "
                           "standard input)",
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
                               "don't write ancestral sequences. IMPORTANT: "
                               "Must be used in conjunction with --refGenome"
                               " to set a non-ancestral genome as the reference"
                               " because the default reference is the root.", 
                               false);
  optionsParser->addOptionFlag("onlySequenceNames",
                               "use only sequence names "
                               "for output names.  By default, the UCSC convention of Genome.Sequence "
                               "is used",
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
  optionsParser->addOption("maxBlockLen",
                           "maximum length of MAF block in output",
                           MafBlock::defaultMaxLength);
  optionsParser->addOptionFlag("global", "output all columns in alignment, "
                               "ignoring refGenome, refSequence, etc. flags",
                               false);
  optionsParser->addOptionFlag("printTree", "print a gene tree for every block",
                               false);
  optionsParser->addOptionFlag("onlyOrthologs", "make only orthologs to the "
                               "reference appear in the MAF blocks", false);

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
  bool global;
  bool printTree;
  bool onlyOrthologs;
  hal_index_t maxBlockLen;
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
    ucscNames = !optionsParser->getFlag("onlySequenceNames");
    unique = optionsParser->getFlag("unique");
    append = optionsParser->getFlag("append");
    global = optionsParser->getFlag("global");
    printTree = optionsParser->getFlag("printTree");
    maxBlockLen = optionsParser->getOption<hal_index_t>("maxBlockLen");
    onlyOrthologs = optionsParser->getFlag("onlyOrthologs");

    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
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

    if (noAncestors == true && refGenome->getNumChildren() != 0 && !global)
    {
      throw hal_exception(string("Since the reference genome to be used for the"
                                 " MAF is ancestral (") + refGenome->getName() +
                          "), the --noAncestors option is invalid.  The "
                          "--refGenome option can be used to specify a "
                          "different reference.");
    }
    
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
    ofstream mafFileStream;
    if (mafPath != "stdout")
    {
      mafFileStream.open(mafPath.c_str(), openFlags);
      if (!mafFileStream)
      {
        throw hal_exception("Error opening " + mafPath);
      }
    }
    ostream& mafStream = mafPath != "stdout" ? mafFileStream : cout;

    MafExport mafExport;
    mafExport.setMaxRefGap(maxRefGap);
    mafExport.setNoDupes(noDupes);
    mafExport.setNoAncestors(noAncestors);
    mafExport.setUcscNames(ucscNames);
    mafExport.setUnique(unique);
    mafExport.setAppend(append);
    mafExport.setMaxBlockLength(maxBlockLen);
    mafExport.setPrintTree(printTree);
    mafExport.setOnlyOrthologs(onlyOrthologs);

    ifstream refTargetsStream;
    if (refTargetsPath != "\"\"")
    {
      ifstream bedFileStream;
      if (refTargetsPath != "stdin")
      {
        bedFileStream.open(refTargetsPath.c_str());
        if (!refTargetsStream)
        {
          throw hal_exception("Error opening " + refTargetsPath);
        }
      }
      istream& bedStream = refTargetsPath != "stdin" ? bedFileStream : cin;
      MafBed mafBed(mafStream, alignment, refGenome, refSequence, start,
                    length, targetSet, mafExport);
      mafBed.scan(&bedStream);
    }
    else
    {
      if (start == 0 && length == 0 && ref->getSequenceLength() == 0)
      {
        string refSeqName = 
           refSequence != NULL ? refSequence->getName() : refGenome->getName();
        cerr << "hal2maf: Warning reference sequence " << refSeqName
             << " has zero length.  MAF output will be empty" << endl;
      }
      else if(global)
      {
        mafExport.convertEntireAlignment(mafStream, alignment);
      }
      else
      {
        mafExport.convertSegmentedSequence(mafStream, alignment, ref, 
                                           start, length, targetSet);
      }
    }
    if (mafPath != "stdout")
    {
      // dont want to leave a size 0 file when there's not ouput because
      // it can make some scripts (ie that process a maf for each contig)
      // obnoxious (presently the case for halPhlyoPTrain which uses 
      // hal2mafMP --splitBySequence). 
      if (mafFileStream.tellp() == (streampos)0)
      {
        std::remove(mafPath.c_str());
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
