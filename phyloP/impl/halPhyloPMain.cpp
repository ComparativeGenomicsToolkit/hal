/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "halPhyloP.h"
#undef min
using namespace std;
using namespace hal;

/** This is basically copied from halAlignability.cpp but returns
 * the phyloP score rather than the count of alignments.
 *
 * Coordinates are always genome-relative by default (as opposed to 
 * sequence-relative).  The one exception is all methods within the
 * Sequence interface. 
 *
 * By default, all bases in the reference genome are scanned.  And all
 * other genomes are considered.  The --refSequence, --start, and 
 * --length options can limit the query to a subrange.  Note that unless
 * --refSequence is specified, --start is genome-relative (based on 
 * all sequences being concatenated together).
 *
 * Genomes to include are determined by genomes in the neutral model 
 * (.mod) file- this should be in the format outputted by phyloFit.
 *
 * There are a few options for dealing with duplications. The default
 * is dupMask=soft, dupType=ambiguous. This replaces any species base 
 * with an N if duplications cause uncertainty in the base for a 
 * particular alignment column. If dupType=all, then all duplications
 * are masked regardless of whether they cause uncertainty. If
 * dupMask=hard, then the entire column is masked (p-value returned is 1.0)
 */

/** Print the phyloP wiggle for a subrange of a given sequence to
 * the output stream. */
static void printSequence(ostream& outStream, halPhyloP *phyloP, 
                          const Sequence* sequence, 
                          hal_size_t start, hal_size_t length, hal_size_t step);

/** If given genome-relative coordinates, map them to a series of 
 * sequence subranges */
static void printGenome(ostream& outStream, halPhyloP *phyloP,
                        const Genome* genome, const Sequence* sequence,
                        hal_size_t start, hal_size_t length, hal_size_t step);

static const hal_size_t StringBufferSize = 1024;

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
  hal_size_t start;
  hal_size_t length;
  hal_size_t step;
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
    

    halPhyloP phyloPObj = halPhyloP(alignment, modPath);
    halPhyloP *phyloP = &phyloPObj;
    phyloP->setDupMask(dupMask);
    phyloP->setDupType(dupType);

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
    
    printGenome(outStream, phyloP, refGenome, refSequence, start, length, 
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

/** Given a Sequence (chromosome) and a (sequence-relative) coordinate
 * range, print the phyloP wiggle with respect to the genomes
 * in the target set */
void printSequence(ostream& outStream, halPhyloP *phyloP, 
		   const Sequence* sequence, 
                   hal_size_t start, hal_size_t length, hal_size_t step)
{
  hal_size_t seqLen = sequence->getSequenceLength();
  if (seqLen == 0)
  {
    return;
  }
  /** If the length is 0, we do from the start position until the end
   * of the sequence */
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

  /** The ColumnIterator is fundamental structure used in this example to
   * traverse the alignment.  It essientially generates the multiple alignment
   * on the fly according to the given reference (in this case the target
   * sequence).  Since this is the sequence interface, the positions
   * are sequence relative.  Note that we must specify the last position
   * in advance when we get the iterator.  This will limit it following
   * duplications out of the desired range while we are iterating. */
  hal_size_t pos = start;
  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(&phyloP->_targetSet,
                                                             0, pos,
                                                             last - 1,
                                                             true);
  // note wig coordinates are 1-based for some reason so we shift to right
  outStream << "fixedStep chrom=" << sequenceName << " start=" << start + 1
            << " step=" << step << "\n";

  // set float precision to 3 places to be consistent with non-hal phyloP
  outStream.setf(ios::fixed, ios::floatfield);
  outStream.precision(3);
  
  /** Since the column iterator stores coordinates in Genome coordinates
   * internally, we have to switch back to genome coordinates.  */
  // convert to genome coordinates
  pos += sequence->getStartPosition();
  last += sequence->getStartPosition();
  while (pos <= last)
  {
    /** ColumnIterator::ColumnMap maps a Sequence to a list of bases
     * the bases in the map form the alignment column.  Some sequences
     * in the map can have no bases (for efficiency reasons) */ 
    const ColumnIterator::ColumnMap* cmap = colIt->getColumnMap();
    double pval = phyloP->pval(cmap);

    outStream << pval << '\n';
    
    /** lastColumn checks if we are at the last column (inclusive)
     * in range.  So we need to check at end of iteration instead
     * of beginning (which would be more convenient).  Need to 
     * merge global fix from other branch */
    if (colIt->lastColumn() == true)
    {
      break;
    }

    pos += step;    
    if (step == 1)
    {
      /** Move the iterator one position to the right */
      colIt->toRight();
      
      /** This is some tuning code that will probably be hidden from 
       * the interface at some point.  It is a good idea to use for now
       * though */
      // erase empty entries from the column.  helps when there are 
      // millions of sequences (ie from fastas with lots of scaffolds)
      if (pos % 1000 == 0)
      {
        colIt->defragment();
      }
    }
    else
    {
      /** Reset the iterator to a non-contiguous position */
      colIt->toSite(pos, last);
    }
  }
}

/** Map a range of genome-level coordinates to potentially multiple sequence
 * ranges.  For example, if a genome contains two chromosomes ChrA and ChrB,
 * both of which are of length 500, then the genome-coordinates would be
 * [0,499] for ChrA and [500,999] for ChrB. All aspects of the HAL API
 * use these global coordinates (chromosomes concatenated together) except
 * for the hal::Sequence interface.  We can convert between the two by 
 * adding or subtracting the sequence start position (in the example it woudl
 * be 0 for ChrA and 500 for ChrB) */
void printGenome(ostream& outStream, halPhyloP *phyloP, 
                 const Genome* genome, const Sequence* sequence,
                 hal_size_t start, hal_size_t length, hal_size_t step)
{
  if (sequence != NULL)
  {
    printSequence(outStream, phyloP, sequence, start, length, step);
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

        printSequence(outStream, phyloP, sequence, readStart, readLen, step);
        runningLength += readLen;
      }
    }
  }
}

