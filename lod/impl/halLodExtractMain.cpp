/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include "halLodNode.h"
#include "halLodEdge.h"
#include "halLodAdjTable.h"
#include "halLodGraph.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("step", "step size for interpolation");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (root if empty)", 
                           "\"\"");
                           
  optionsParser->setDescription("Extract LOD");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();
  string halPath;
  string refGenomeName;
  hal_size_t step;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    step = optionsParser->getArgument<hal_size_t>("step");
    refGenomeName = optionsParser->getOption<string>("refGenome");
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
    if (refGenomeName == "\"\"")
    {
      refGenomeName = alignment->getRootName();
    }
    const Genome* parent = alignment->openGenome(refGenomeName);
    if (parent == NULL)
    {
      throw hal_exception(string("Genome ") + refGenomeName + " not found");
    }

    vector<const Genome*> children;
    for (hal_size_t child = 0; child < parent->getNumChildren(); ++child)
    {
      children.push_back(parent->getChild(child));
    }

    LodGraph lodGraph;
    lodGraph.build(alignment, parent, children, step);
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
