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
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("inHalPath", "Input hal file");
  optionsParser->addArgument("outHalPath", "output hal file");
  optionsParser->addArgument("step", "Step size for interpolation");
  optionsParser->addOption("root", 
                           "Name of root genome of tree to extract (root if "
                           "empty)", "\"\"");
  optionsParser->addOption("outTree", 
                           "Phylogenetic tree of output HAL file.  Must "
                           "contain only genomes from the input HAL file. "
                           "(input\'s tree used if empty)", "\"\"");
  optionsParser->addOptionFlag("keepSequences",
                               "Write the sequence strings to the output "
                               "file.", false);
  optionsParser->setDescription("Generate a new HAL file at a coarser "
                                "Level of Detail (LOD) by interpolation.");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string inHalPath;
  string outHalPath;
  string rootName;
  string outTree;
  hal_size_t step;
  bool keepSequences;
  try
  {
    optionsParser->parseOptions(argc, argv);
    inHalPath = optionsParser->getArgument<string>("inHalPath");
    outHalPath = optionsParser->getArgument<string>("outHalPath");
    rootName = optionsParser->getOption<string>("root");
    outTree = optionsParser->getOption<string>("outTree");
    step = optionsParser->getArgument<hal_size_t>("step");
    keepSequences = optionsParser->getFlag("keepSequences");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
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
                                           step, outTree, rootName,
                                           keepSequences);
  }
/*  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
*/
  catch(float temp) {}
  return 0;
}
