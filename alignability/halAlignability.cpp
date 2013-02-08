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

using namespace std;
using namespace hal;


static void printSequence(ostream& outStream, const Sequence* sequence, 
                          const set<const Genome*>& targetSet,
                          hal_size_t start, hal_size_t length, hal_size_t step);
static void printGenome(ostream& outStream,
                        const Genome* genome, const Sequence* sequence,
                        const set<const Genome*>& targetSet,
                        hal_size_t start, hal_size_t length, hal_size_t step);

static const hal_size_t StringBufferSize = 1024;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(false);
  optionsParser->addArgument("halPath", "input hal file");
  optionsParser->addArgument("refGenome", "reference genome to scan");
  optionsParser->addOption("outWiggle", "output wig file (stdout if none)",
                           "stdout");
  optionsParser->addOption("refSequence", "sequence name to export ("
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
  optionsParser->addOption("rootGenome", 
                           "name of root genome (none if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOption("step", "step size", 1);
  optionsParser->setDescription("Make alignability wiggle plot for a genome.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string wigPath;
  string refGenomeName;
  string rootGenomeName;
  string targetGenomes;
  string refSequenceName;
  hal_size_t start;
  hal_size_t length;
  hal_size_t step;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halPath");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    wigPath = optionsParser->getOption<string>("outWiggle");
    refSequenceName = optionsParser->getOption<string>("refSequence");
    start = optionsParser->getOption<hal_size_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    rootGenomeName = optionsParser->getOption<string>("rootGenome");
    targetGenomes = optionsParser->getOption<string>("targetGenomes");
    step = optionsParser->getOption<hal_size_t>("step");

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
      throw hal_exception("input hal alignmenet is empty");
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

    ofstream ofile;
    ostream& outStream = wigPath == "stdout" ? cout : ofile;
    if (wigPath != "stdout")
    {
      ofile.open(wigPath.c_str());
      if (!ofile)
      {
        throw hal_exception(string("Error opening output file ") + 
                            wigPath);
      }
    }
    
    printGenome(outStream, refGenome, refSequence, targetSet, start, length, 
                step);
    
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

void printSequence(ostream& outStream, const Sequence* sequence, 
                   const set<const Genome*>& targetSet,
                   hal_size_t start, hal_size_t length, hal_size_t step)
{
  hal_size_t seqLen = sequence->getSequenceLength();
  if (length == 0)
  {
    length = seqLen - start;
  }
  hal_size_t last = start + length;
  if (last > seqLen)
  {
    stringstream ss;
    ss << "Specified range [" << start << "," << length << "] is"
       << "out of range for sequence " << sequence->getName() 
       << ", which has length " << seqLen;
    throw (hal_exception(ss.str()));
  }

  const Genome* genome = sequence->getGenome();
  string sequenceName = sequence->getName();
  string genomeName = genome->getName();
  if (sequenceName.find(genomeName + '.') == 0)
  {
    sequenceName = sequenceName.substr(genomeName.length() + 1);
  }

  hal_size_t pos = start;
  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(&targetSet,
                                                             0, pos,
                                                             last - 1,
                                                             true);
  outStream << "fixedStep\tchrom=" << sequenceName << "\tstart=" << start
            << "\tstep=" << step << "\n";
  
  // convert to genome coordinates
  pos += sequence->getStartPosition();
  last += sequence->getStartPosition();
  while (colIt->lastColumn() == false && pos <= last)
  {
    hal_size_t count = 0;
    const ColumnIterator::ColumnMap* cmap = colIt->getColumnMap();
    for (ColumnIterator::ColumnMap::const_iterator i = cmap->begin();
         i != cmap->end(); ++i)
    {
      if (i->first->getGenome() != genome && !i->second->empty())
      {
        ++count; 
      }
    }

    outStream << count << '\n';

    pos += step;    
    if (step == 1)
    {
      colIt->toRight();
      // erase empty entries from the column.  helps when there are 
      // millions of sequences (ie from fastas with lots of scaffolds)
      if (pos % 1000 == 0)
      {
        colIt->defragment();
      }
    }
    else
    {
      colIt->toSite(pos, last);
    }
  }
}

void printGenome(ostream& outStream,
                 const Genome* genome, const Sequence* sequence,
                 const set<const Genome*>& targetSet,
                 hal_size_t start, hal_size_t length, hal_size_t step)
{
  if (sequence != NULL)
  {
    printSequence(outStream, sequence, targetSet, start, length, step);
  }
  else
  {
    if (start + length > genome->getSequenceLength())
    {
      stringstream ss;
      ss << "Specified range [" << start << "," << length << "] is"
         << "out of range for genome " << genome->getName() 
         << ", which has length " << genome->getSequenceLength();
      throw (hal_exception(ss.str()));
    }
    if (length == 0)
    {
      length = genome->getSequenceLength() - start;
    }

    SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
    SequenceIteratorConstPtr seqEndIt = genome->getSequenceEndIterator();
    hal_size_t runningLength = 0;
    for (; seqIt != seqEndIt; seqIt->toNext())
    {
      const Sequence* sequence = seqIt->getSequence();
      hal_size_t seqLen = sequence->getSequenceLength();
      hal_size_t seqStart = (hal_size_t)sequence->getStartPosition();

      if (start + length >= seqStart && 
          start < seqStart + seqLen &&
          runningLength < length)
      {
        hal_size_t readStart = seqStart >= start ? 0 : seqStart - start;
        hal_size_t readLen = std::min(seqLen - start, length - runningLength);

        printSequence(outStream, sequence, targetSet, readStart, readLen, step);
        runningLength += readLen;
      }
    }
  }
}

