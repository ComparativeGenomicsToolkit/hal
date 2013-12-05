/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu) and 
 * Melissa Jane Hubisz (Cornell University)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "halPhyloP.h"
#include "halPhyloPBed.h"

#undef min
using namespace std;
using namespace hal;

/** If given genome-relative coordinates, map them to a series of 
 * sequence subranges */
static void printGenome(PhyloP *phyloP,
                        const Genome* genome, const Sequence* sequence,
                        hal_size_t start, hal_size_t length, hal_size_t step);


static CLParserPtr initParser()
{
  /** It is convenient to use the HAL command line parser for the command
   * line because it automatically adds some comman options.  Using the 
   * parser is by no means required however */
  CLParserPtr optionsParser = hdf5CLParserInstance(false);
  optionsParser->addArgument("halPath", "input hal file");
  optionsParser->addArgument("refGenome", "reference genome to scan");
  optionsParser->addArgument("modPath", "input neutral model file");
  optionsParser->addArgument("outWiggle", "output wig file (or \"stdout\""
                             " to pipe to standard output)");
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
  optionsParser->addOption("dupType",
                           "Which duplications to mask according to dupMask "
                           "option. Choices are: "
                           "\"all\": Any duplicated region; or "
                           "\"ambiguous\": Regions within duplications where "
                           "alignments from the same species do not contain"
                           " the same base.",
                           "ambiguous");
  optionsParser->addOption("dupMask",
                           "What to do with duplicated regions. Choices are: "
                           "\"hard\": mask entire alignment column if any "
                           "duplications occur; or "
                           "\"soft\": mask species where duplications occur.",
                           "soft");
  optionsParser->addOption("step", "step size", 1);
  optionsParser->addOption("refBed", "Bed file with coordinates to "
                           "annotate in the reference genome (or \"stdin\") "
                           "to stream from standard input", "\"\"");
  optionsParser->addOption("subtree", "Partition the tree into the subtree "
			   "beneath the node given (including the branch "
			   "leading up to this node), and test for "
			   "conservation/acceleration in this subtree "
			   "relative to the rest of the tree", "\"\"");
  optionsParser->addOption("prec", "Number of decimal places in wig output", 3);
  
  optionsParser->setDescription("Make PhyloP wiggle plot for a genome.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string modPath;
  string halPath;
  string wigPath;
  string refGenomeName;
  string refSequenceName;
  string dupType;
  string dupMask;
  string subtree;
  hal_size_t start;
  hal_size_t length;
  hal_size_t step;
  string refBedPath;
  hal_size_t prec;
  try
  {
    optionsParser->parseOptions(argc, argv);
    modPath = optionsParser->getArgument<string>("modPath");
    halPath = optionsParser->getArgument<string>("halPath");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    wigPath = optionsParser->getArgument<string>("outWiggle");
    refSequenceName = optionsParser->getOption<string>("refSequence");
    start = optionsParser->getOption<hal_size_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    step = optionsParser->getOption<hal_size_t>("step");
    dupType = optionsParser->getOption<string>("dupType");
    dupMask = optionsParser->getOption<string>("dupMask");
    subtree = optionsParser->getOption<string>("subtree");
    std::transform(dupMask.begin(), dupMask.end(), dupMask.begin(), ::tolower);
    refBedPath = optionsParser->getOption<string>("refBed");
    prec = optionsParser->getOption<hal_size_t>("prec");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  try
  {
    /** Everything begins with the alignment object, which is created
     * via a path to a .hal file.  Options don't necessarily need to
     * come from the optionsParser -- see other interfaces in 
     * hal/api/inc/halAlignmentInstance.h */
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, 
                                                           optionsParser);

    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("input hal alignmenet is empty");
    }

    /** Open the reference genome */
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
    
    /** If a sequence was spefied we look for it in the reference genome */
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
    
    // set the precision of floating point output
    outStream.setf(ios::fixed, ios::floatfield);
    outStream.precision(prec);
    
    PhyloP phyloP;
    phyloP.init(alignment, modPath, &outStream, dupMask == "soft" , dupType,
		"CONACC", subtree);

    ifstream refBedStream;
    if (refBedPath != "\"\"")
    {
      ifstream bedFileStream;
      if (refBedPath != "stdin")
      {
        bedFileStream.open(refBedPath.c_str());
        if (!refBedStream)
        {
          throw hal_exception("Error opening " + refBedPath);
        }
      }
      istream& bedStream = refBedPath != "stdin" ? bedFileStream : cin;
      PhyloPBed phyloPBed(alignment, refGenome, refSequence, 
                          start, length, step, phyloP, outStream);
      phyloPBed.scan(&bedStream);
    }
    else
    {
      printGenome(&phyloP, refGenome, refSequence, start, length, step);
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


/** Map a range of genome-level coordinates to potentially multiple sequence
 * ranges.  For example, if a genome contains two chromosomes ChrA and ChrB,
 * both of which are of length 500, then the genome-coordinates would be
 * [0,499] for ChrA and [500,999] for ChrB. All aspects of the HAL API
 * use these global coordinates (chromosomes concatenated together) except
 * for the hal::Sequence interface.  We can convert between the two by 
 * adding or subtracting the sequence start position (in the example it woudl
 * be 0 for ChrA and 500 for ChrB) */
void printGenome(PhyloP *phyloP, 
                 const Genome* genome, const Sequence* sequence,
                 hal_size_t start, hal_size_t length, hal_size_t step)
{
  if (sequence != NULL)
  {
    phyloP->processSequence(sequence, start, length, step);
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
        hal_size_t readStart = seqStart >= start ? 0 : start - seqStart;
        hal_size_t readLen = min(seqLen - readStart, length);
        readLen = min(readLen, length - runningLength);
        phyloP->processSequence(sequence, readStart, readLen, step);
        runningLength += readLen;
      }
    }
  }
}

