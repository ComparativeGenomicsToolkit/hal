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

static void extractAlignedRegions(const Genome* genome,
                                  ostream* bedStream,
                                  bool complement,
                                  bool viewParentCoords);

static void initParser(CLParser& optionsParser) {
  optionsParser.setDescription("Extract aligned regions of genome (with its"
                                " parent) and output to bed file");
  optionsParser.addArgument("halPath", "input hal file");
  optionsParser.addArgument("genome", "name of genome to scan"
                             " (can't be root)");
  optionsParser.addOption("alignedFile", "path to bed file to write to", 
                           "stdout");
  optionsParser.addOptionFlag("complement", "extract the regions of the "
                               "genome that are *unaligned* to the parent. "
                               "ie all intervals that are not returned with"
                               " the default setting.", false);
  optionsParser.addOptionFlag("viewParentCoords", "view the corresponding parent "
                               "coordinates for aligned regions", false);

}

int main(int argc, char** argv)
{
    CLParser optionsParser;
    initParser(optionsParser);

  string halPath;
  string genomeName;
  string outBedPath;
  bool complement;
  bool viewParentCoords;
  try
  {
    optionsParser.parseOptions(argc, argv);
    halPath = optionsParser.getArgument<string>("halPath");
    genomeName = optionsParser.getArgument<string>("genome");
    outBedPath = optionsParser.getOption<string>("alignedFile");
    complement = optionsParser.getFlag("complement");
    viewParentCoords = optionsParser.getFlag("viewParentCoords");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser.printUsage(cerr);
    exit(1);
  }

  try
  {
      AlignmentConstPtr inAlignment(openHalAlignment(halPath, &optionsParser));
    if (inAlignment->getNumGenomes() == 0)
    {
      throw hal_exception("input hal alignment is empty");
    }
    
    const Genome* genome = inAlignment->openGenome(genomeName);
    
    if (genome == NULL)
    {
      throw hal_exception(string("Unable to open genome ") + genomeName +
                          "in alignment.");
    }
    if (genome->getParent() != NULL)
    {
      ostream* bedStream = &cout;
      if (outBedPath != "stdout")
      {
        bedStream = new ofstream(outBedPath.c_str());
      }
      if (!bedStream)
      {
        throw hal_exception(string("Error opening ") + outBedPath + 
                            " for writing");
      }
      extractAlignedRegions(genome, bedStream, complement, viewParentCoords);
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

void extractAlignedRegions(const Genome* genome, ostream* bedStream,
                           bool complement, bool viewParentCoords)
{
  assert(genome && bedStream);
  if (genome->getNumTopSegments() > 0)
  {
    TopSegmentIteratorPtr topSeg = genome->getTopSegmentIterator();
    for (; (not topSeg->atEnd()); topSeg->toRight())
    {
      if ((topSeg->tseg()->hasParent() && complement == false) ||
          (!topSeg->tseg()->hasParent() && complement == true))
      {
        const Sequence* sequence = topSeg->getSequence();
        *bedStream << sequence->getName() << '\t'
                   << topSeg->getStartPosition() - sequence->getStartPosition()
                   << '\t' 
                   << (1 + topSeg->getEndPosition() -
                       sequence->getStartPosition());
        if (!complement && viewParentCoords) {
            BottomSegmentIteratorPtr botSeg = genome->getParent()->getBottomSegmentIterator();
            botSeg->toParent(topSeg);
            if (botSeg->getReversed()) {
                // Ensure the parent segment coordinates always have start < end.
                botSeg->toReverse();
            }
            const Sequence *parentSequence = botSeg->getSequence();
            *bedStream << '\t'
                       << parentSequence->getName() << '\t'
                       << botSeg->getStartPosition() - parentSequence->getStartPosition() << '\t'
                       << 1 + botSeg->getEndPosition() - parentSequence->getStartPosition() << '\t'
                       << (topSeg->tseg()->getParentReversed() ? "-" : "+");
        }
        *bedStream << '\n';
      }
    }
  }
  *bedStream << endl;
}
