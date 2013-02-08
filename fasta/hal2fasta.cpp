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

static void printSequenceLine(ostream& outStream, const Sequence* sequence,
                              hal_size_t start, hal_size_t length, 
                              string& buffer);
static void printSequence(ostream& outStream, const Sequence* sequence, 
                          hal_size_t lineWidth, 
                          hal_size_t start, hal_size_t length);
static void printGenome(ostream& outStream,
                        const Genome* genome, const Sequence* sequence,
                        hal_size_t lineWidth, 
                        hal_size_t start, hal_size_t length);

static const hal_size_t StringBufferSize = 1024;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(false);
  optionsParser->addArgument("inHalPath", "input hal file");
  optionsParser->addArgument("genome", "genome to export");
  optionsParser->addOption("outFaPath", "output fasta file (stdout if none)",
                           "stdout");
  optionsParser->addOption("lineWidth", "Line width for output", 80);
  optionsParser->addOption("sequence", "sequence name to export ("
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
  optionsParser->setDescription("Export single genome from hal database to "
                                "fasta file.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string faPath;
  hal_size_t lineWidth;
  string genomeName;
  string sequenceName;
  hal_size_t start;
  hal_size_t length;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("inHalPath");
    genomeName = optionsParser->getArgument<string>("genome");
    faPath = optionsParser->getOption<string>("outFaPath");
    lineWidth = optionsParser->getOption<hal_size_t>("lineWidth");
    sequenceName = optionsParser->getOption<string>("sequence");
    start = optionsParser->getOption<hal_size_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
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
    
    const Genome* genome = alignment->openGenome(genomeName);
    if (genome == NULL)
    {
      throw hal_exception(string("Genome ") + genomeName + " not found");
    }

    const Sequence* sequence = NULL;
    if (sequenceName != "\"\"")
    {
      sequence = genome->getSequence(sequenceName);
      if (sequence == NULL)
      {
        throw hal_exception(string("Sequence ") + sequenceName + " not found");
      }
    }

    ofstream ofile;
    ostream& outStream = faPath == "stdout" ? cout : ofile;
    if (faPath != "stdout")
    {
      ofile.open(faPath.c_str());
      if (!ofile)
      {
        throw hal_exception(string("Error opening output file ") + 
                            faPath);
      }
    }
    
    printGenome(outStream, genome, sequence, lineWidth, start, length);
    
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

void printSequenceLine(ostream& outStream, const Sequence* sequence,
                       hal_size_t start, hal_size_t length, 
                       string& buffer)
{
  hal_size_t last = start + length;
  hal_size_t readLen;
  for (hal_size_t i = start; i < last; i += StringBufferSize)
  {
    readLen = std::min(StringBufferSize, last - i);
    sequence->getSubString(buffer, i, readLen);
    outStream << buffer;
  }
}

void printSequence(ostream& outStream, const Sequence* sequence, 
                   hal_size_t lineWidth, 
                   hal_size_t start, hal_size_t length)
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
  outStream << '>' << sequence->getName() << '\n';
  hal_size_t readLen;
  string buffer;
  for (hal_size_t i = start; i < last; i += lineWidth)
  {
    readLen = std::min(lineWidth, last - i);
    printSequenceLine(outStream, sequence, i, readLen, buffer);
    outStream << '\n';
  }
}

void printGenome(ostream& outStream,
                 const Genome* genome, const Sequence* sequence,
                 hal_size_t lineWidth, 
                 hal_size_t start, hal_size_t length)
{
  if (sequence != NULL)
  {
    printSequence(outStream, sequence, lineWidth, start, length);
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

        printSequence(outStream, sequence, lineWidth, readStart, readLen);
        runningLength += readLen;
      }
    }
  }
}

