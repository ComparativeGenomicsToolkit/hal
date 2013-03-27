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
  hal_size_t step;
  try
  {
    optionsParser->parseOptions(argc, argv);
    inHalPath = optionsParser->getArgument<string>("inHalPath");
    outHalPath = optionsParser->getArgument<string>("outHalPath");
    rootName = optionsParser->getOption<string>("root");
    step = optionsParser->getArgument<hal_size_t>("step");
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
    if (rootName == "\"\"")
    {
      rootName = inAlignment->getRootName();
    }
    const Genome* parent = inAlignment->openGenome(rootName);
    if (parent == NULL)
    {
      throw hal_exception(string("Genome ") + rootName + " not found");
    }

    vector<const Genome*> children;
    for (hal_size_t child = 0; child < parent->getNumChildren(); ++child)
    {
      children.push_back(parent->getChild(child));
    }

    LodExtract lodExtract;
    lodExtract.createInterpolatedAlignment(inAlignment, outAlignment,
                                           step, "");
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
