/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include "halLodExtract.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inHalPath", "Input hal file");
  optionsParser->addArgument("outHalPath", "output hal file");
  optionsParser->addArgument("scale", "Scale factor for interpolation");
  optionsParser->addOption("root", 
                           "Name of root genome of tree to extract (root if "
                           "empty)", "\"\"");
  optionsParser->addOption("outTree", 
                           "Phylogenetic tree of output HAL file.  Must "
                           "contain only genomes from the input HAL file. "
                           "(input\'s tree used if empty)", "\"\"");
  optionsParser->addOption("probeFrac", 
                           "Fraction of bases in step-interval to sample while "
                           "looking for most aligned column.",
                           0.035);
  optionsParser->addOption("minSeqFrac", 
                           "Minumum sequence length to sample as fraction of "
                           "step size: ie sequences with length <= floor("
                           "minSeqFrac * step) are ignored.",
                           // Note: needs to be manually synched with 
                           // value in halLodInterpolate.py
                           0.5);
  optionsParser->addOptionFlag("keepSequences",
                               "Write the sequence strings to the output "
                               "file.", false);
  optionsParser->addOptionFlag("allSequences", 
                               "Sample all sequences (chromsomes / contigs /"
                               " etc.) no matter how small they are.  "
                               "By default, small sequences may be skipped if "
                               "they fall within the step size.", false);
  optionsParser->setDescription("Generate a new HAL file at a coarser "
                                "Level of Detail (LOD) by interpolation. "
                                "The scale parameter is used to estimate "
                                "the interpolation step-size so that the "
                                "output has \"scale\" fewer blocks than the"
                                " input.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string inHalPath;
  string outHalPath;
  string rootName;
  string outTree;
  double scale;
  bool keepSequences;
  bool allSequences;
  double probeFrac;
  double minSeqFrac;
  try
  {
    optionsParser->parseOptions(argc, argv);
    inHalPath = optionsParser->getArgument<string>("inHalPath");
    outHalPath = optionsParser->getArgument<string>("outHalPath");
    rootName = optionsParser->getOption<string>("root");
    outTree = optionsParser->getOption<string>("outTree");
    scale = optionsParser->getArgument<double>("scale");
    keepSequences = optionsParser->getFlag("keepSequences");
    allSequences = optionsParser->getFlag("allSequences");
    probeFrac = optionsParser->getOption<double>("probeFrac");
    minSeqFrac = optionsParser->getOption<double>("minSeqFrac");
    if (allSequences == true)
    {
      minSeqFrac = 0.;
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    return 1;
  }
  try
  {
    AlignmentConstPtr inAlignment = openHalAlignmentReadOnly(inHalPath, 
                                                             optionsParser);
    if (inAlignment->getNumGenomes() == 0)
    {
      throw hal_exception("Input hal alignment is empty");
    }

    AlignmentPtr outAlignment = hdf5AlignmentInstance();
    outAlignment->setOptionsFromParser(optionsParser);
    outAlignment->createNew(outHalPath);
    
    if (outAlignment->getNumGenomes() != 0)
    {
      throw hal_exception("Output hal Alignmnent cannot be initialized");
    }
    if (rootName != "\"\"" && inAlignment->openGenome(rootName) == NULL)
    {
      throw hal_exception(string("Genome ") + rootName + " not found");
    }
    if (rootName == "\"\"")
    {
      rootName = "";
    }
    if (outTree == "\"\"")
    {
      outTree = "";
    }

    LodExtract lodExtract;
    lodExtract.createInterpolatedAlignment(inAlignment, outAlignment,
                                           scale, outTree, rootName,
                                           keepSequences, allSequences,
                                           probeFrac, minSeqFrac);
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
