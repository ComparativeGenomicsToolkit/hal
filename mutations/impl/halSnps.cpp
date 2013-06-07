/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "hal.h"

using namespace std;
using namespace hal;

static pair<hal_size_t, hal_size_t>
countSnps(const Genome* refGenome, const Genome* queryGenome, 
                 hal_index_t start, hal_size_t length, bool doDupes, 
                 ostream& refBedStream);

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("refGenome", 
                             "name of reference genome.");
  optionsParser->addArgument("queryGenome", 
                             "name of query genome.");
  optionsParser->addOption("bed",
                           "path of bed file to write snps to in reference "
                           "genome coordinates",
                           "\"\"");
  optionsParser->addOption("noDupes",
                           "do not consider paralogies while mapping",
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
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
  optionsParser->setDescription("Count snps between two genomes.  Outputs"
                                " totalSnps totalAlignedPairs");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string refGenomeName;
  string queryGenomeName;
  string bedPath;
  bool noDupes;
  string refSequenceName;
  hal_index_t start;
  hal_size_t length;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    queryGenomeName = optionsParser->getArgument<string>("queryGenome");
    bedPath = optionsParser->getOption<string>("bed");
    noDupes = optionsParser->getOption<bool>("noDupes");
    refSequenceName = optionsParser->getOption<string>("refSequence");
    start = optionsParser->getOption<hal_index_t>("start");
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
      throw hal_exception("hal alignmenet is empty");
    }
    
    const Genome* queryGenome = NULL;
    queryGenome = alignment->openGenome(queryGenomeName);
    if (queryGenome == NULL)
    {
      throw hal_exception(string("Query genome, ") + queryGenomeName + 
                          ", not found in alignment");
    }

    const Genome* refGenome = NULL;
    refGenome = alignment->openGenome(refGenomeName);
    if (refGenome == NULL)
    {
      throw hal_exception(string("Reference genome, ") + refGenomeName + 
                          ", not found in alignment");
    }

    if (start + length >= refGenome->getSequenceLength())
    {
      throw hal_exception(string("Invalid range for ") + refGenomeName);
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
    if (bedPath != "\"\"")
    {
      refBedStream.open(bedPath.c_str());
      if (!refBedStream)
      {
        throw hal_exception("Error opening " + bedPath);
      }
    }

    pair<hal_size_t, hal_size_t> res = countSnps(refGenome, queryGenome, 
                                                 start, length, !noDupes,
                                                 refBedStream);
    cout << res.first << " " << res.second << endl;
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

pair<hal_size_t, hal_size_t>
countSnps(const Genome* refGenome, const Genome* queryGenome, 
          hal_index_t start, hal_size_t length, bool doDupes, 
          ostream& refBedStream)
{
  SegmentIteratorConstPtr refIt;
  hal_index_t lastIndex = 0;
  hal_index_t lastPos = start + length - 1;
  if (refGenome->getNumTopSegments() > 0)
  {
    refIt = refGenome->getTopSegmentIterator();
    lastIndex = (hal_index_t)refGenome->getNumTopSegments();
  }
  else if (refGenome->getNumBottomSegments() > 0)

  {
    refIt = refGenome->getBottomSegmentIterator();
    lastIndex = (hal_index_t)refGenome->getNumBottomSegments();
  }
  if (refIt.get() == NULL)
  {
    return pair<hal_size_t, hal_size_t>(0,0);
  }
  // move segment to just the start position
  refIt->toSite(start, true);
  // unslice right of start position
  refIt->slice(refIt->getStartOffset(), 0);
  // but slice if we need to cut by length
  if (length < refIt->getLength())
  {
    refIt->slice(refIt->getStartOffset(), refIt->getLength() - length);
  }
    
  set<const Genome*> inputSet;
  inputSet.insert(refGenome);
  inputSet.insert(queryGenome);
  set<const Genome*> spanningTree;
  getGenomesInSpanningTree(inputSet, spanningTree);
  string refString;
  string queryString;
  hal_size_t subCount = 0;
  hal_size_t colCount = 0;

  set<hal::MappedSegmentConstPtr> mappedSet;
    
  while (refIt->getArrayIndex() < lastIndex &&
         refIt->getStartPosition() <= lastPos)  
  {
    refIt->getMappedSegments(mappedSet, queryGenome, &spanningTree,
                              doDupes, 0);
    refIt->toRight(lastPos);
  }
  
  for (set<MappedSegmentConstPtr>::const_iterator i = mappedSet.begin();
       i != mappedSet.end(); ++i)
  {
    SlicedSegmentConstPtr refSeg = (*i)->getSource();
    (*i)->getString(queryString);
    refSeg->getString(refString);
    const Sequence* refSequence = refSeg->getSequence();
    hal_index_t refSeqStart = refSequence->getStartPosition();
    colCount += refString.length();
    for (size_t j = 0; j < refString.length(); ++j)
    {
      if (isSubstitution(refString[j], queryString[j]))
      {
        ++subCount;
        if (refBedStream)
        {
          refBedStream << refSequence->getName() << '\t' 
                       << refSeg->getStartPosition() + j - refSeqStart << '\t'
                       << 1 + refSeg->getStartPosition() + j - refSeqStart
                       << '\t' << refString[j] << queryString[j] << "\t0\t"
                       << ((*i)->getReversed() ? '-' : '+') << '\n';
        }
      } 
    }
  }
  return pair<hal_size_t, hal_size_t>(subCount, colCount);
}
