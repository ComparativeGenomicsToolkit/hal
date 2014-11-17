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

/** This is a tool that counts the number of other genomes each base in
 * a query region is aligned to.
 *
 * Coordinates are always genome-relative by default (as opposed to 
 * sequence-relative).  The one exception is all methods within the
 * Sequence interface. 
 *
 * By default, all bases in the referecen genome are scanned.  And all
 * other genomes are considered.  The --refSequence, --start, and 
 * --length options can limit the query to a subrange.  Note that unless
 * --refSequence is specified, --start is genome-relative (based on 
 * all sequences being concatenated together).
 *
 * Other genomes to query (default all) are controlled by --rootGenome
 * (name of highest ancestor to consider) and/or --targetGenomes
 * (a list of genomes to consider).  
 *
 * So if a base in the reference genome is aligned to a base in a genome
 * that is not under root or in the target list, it will not count to the
 * alignment depth.
 */

/** Print the alignment depth wiggle for a subrange of a given sequence to
 * the output stream. */
static void printSequence(ostream& outStream, const Sequence* sequence, 
                          const set<const Genome*>& targetSet,
                          hal_size_t start, hal_size_t length, hal_size_t step,
                          bool countDupes, bool noAncestors);

/** If given genome-relative coordinates, map them to a series of 
 * sequence subranges */
static void printGenome(ostream& outStream,
                        const Genome* genome, const Sequence* sequence,
                        const set<const Genome*>& targetSet,
                        hal_size_t start, hal_size_t length, hal_size_t step,
                        bool countDupes, bool noAncestors);

static const hal_size_t StringBufferSize = 1024;

static CLParserPtr initParser()
{
  /** It is convenient to use the HAL command line parser for the command
   * line because it automatically adds some comman options.  Using the 
   * parser is by no means required however */
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
  optionsParser->addOptionFlag("countDupes",
                               "count each other *position* each base aligns "
                               "to, rather than the number of unique genomes, "
                               "including paralogies so a genome can be "
                               "counted  multiple times.  This will give the "
                               "height of the MAF column created with hal2maf.",
                               false);
  optionsParser->addOptionFlag("noAncestors", 
                               "do not count ancestral genomes.", false);
  optionsParser->setDescription("Make alignment depth wiggle plot for a genome. "
                                "By default, this is a count of the number of "
                                "other unique genomes each base aligns to, "
                                "including ancestral genomes.");
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
  bool countDupes;
  bool noAncestors;
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
    countDupes = optionsParser->getFlag("countDupes");
    noAncestors = optionsParser->getFlag("noAncestors");

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
    
    /** Alignments are composed of sets of Genomes.  Each genome is a set
     * of Sequences (chromosomes).  They are accessed by their names.  
     * here we map the root and targetSet parameters (if specifeid) to 
     * a sset of readonly Genome pointers */
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
    const SegmentedSequence* ref = refGenome;
    
    /** If a sequence was spefied we look for it in the reference genome */
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

    if (refGenome->getNumChildren() != 0 && noAncestors == true)
    {
      throw hal_exception(string("--noAncestors cannot be used when reference "
                                 "genome (") + refGenome->getName() + 
                          string(") is ancetral"));
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
                step, countDupes, noAncestors);
    
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
 * range, print the alignmability wiggle with respect to the genomes
 * in the target set */
void printSequence(ostream& outStream, const Sequence* sequence, 
                   const set<const Genome*>& targetSet,
                   hal_size_t start, hal_size_t length, hal_size_t step,
                   bool countDupes, bool noAncestors)
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
  ColumnIteratorConstPtr colIt = sequence->getColumnIterator(&targetSet,
                                                             0, pos,
                                                             last - 1,
                                                             false,
                                                             noAncestors);
  // note wig coordinates are 1-based for some reason so we shift to right
  outStream << "fixedStep chrom=" << sequenceName << " start=" << start + 1
            << " step=" << step << "\n";
  
  /** Since the column iterator stores coordinates in Genome coordinates
   * internally, we have to switch back to genome coordinates.  */
  // convert to genome coordinates
  pos += sequence->getStartPosition();
  last += sequence->getStartPosition();
  // keep track of unique genomes
  set<const Genome*> genomeSet;
  while (pos <= last)
  {
    genomeSet.clear();
    hal_size_t count = 0;
    /** ColumnIterator::ColumnMap maps a Sequence to a list of bases
     * the bases in the map form the alignment column.  Some sequences
     * in the map can have no bases (for efficiency reasons) */ 
    const ColumnIterator::ColumnMap* cmap = colIt->getColumnMap();

    /** For every sequence in the map */
    for (ColumnIterator::ColumnMap::const_iterator i = cmap->begin();
         i != cmap->end(); ++i)
    {
      if (countDupes == true)
      {
        // countDupes enabled: we just count everything
        count += i->second->size();
      }
      else if (!i->second->empty())
      {
        // just counting unique genomes: add it if there's at least one base
        genomeSet.insert(i->first->getGenome());
      }
    }
    if (countDupes == false) 
    {
      count = genomeSet.size();
    }
    // don't want to include reference base in output
    --count;

    outStream << count << '\n';
    
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
void printGenome(ostream& outStream,
                 const Genome* genome, const Sequence* sequence,
                 const set<const Genome*>& targetSet,
                 hal_size_t start, hal_size_t length, hal_size_t step,
                 bool countDupes, bool noAncestors)
{
  if (sequence != NULL)
  {
    printSequence(outStream, sequence, targetSet, start, length, step, countDupes,
                  noAncestors);
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
        printSequence(outStream, sequence, targetSet, readStart, readLen, step,
                      countDupes, noAncestors);
        runningLength += readLen;
      }
    }
  }
}

