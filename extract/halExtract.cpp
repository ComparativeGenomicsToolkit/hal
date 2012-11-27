/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include "hal.h"

using namespace std;
using namespace hal;

static void getDimensions(const Genome* genome,
                          vector<Sequence::Info>& dimensions);

static void extract(const AlignmentConstPtr inAlignment,
                    AlignmentPtr outAlignment, const std::string& rootName);

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("inHalPath", "input hal file");
  optionsParser->addArgument("outHalPath", "output hal file");
  optionsParser->addOption("root", "root of subtree to extract", "\"\"");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string inHalPath;
  string outHalPath;
  string rootName;
  try
  {
    optionsParser->parseOptions(argc, argv);
    inHalPath = optionsParser->getArgument<string>("inHalPath");
    outHalPath = optionsParser->getArgument<string>("outHalPath");
    rootName = optionsParser->getOption<string>("root");
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
      throw hal_exception("input hal alignmenet is empty");
    }

    AlignmentPtr outAlignment = openHalAlignment(outHalPath, 
                                                 optionsParser);
    outAlignment->createNew(outHalPath);
    if (outAlignment->getNumGenomes() != 0)
    {
      throw hal_exception("output hal alignmenet cannot be initialized");
    }
    
    if (rootName != "\"\"" && inAlignment->getNumGenomes() != 0)
    {
      rootName = inAlignment->getRootName();
    }

    extract(inAlignment, outAlignment, rootName);
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

void getDimensions(const Genome* genome, vector<Sequence::Info>& dimensions)
{

}

void extract(const AlignmentConstPtr inAlignment,
             AlignmentPtr outAlignment, const std::string& rootName)
{
  const Genome* genome = inAlignment->openGenome(rootName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome not found: ") + rootName);
  }
  Genome* newGenome = NULL;
  if (outAlignment->getNumGenomes() == 0 || genome->getParent() == NULL)
  {
    newGenome = outAlignment->addRootGenome(rootName);
  }
  else
  {
    const Genome* parent = genome->getParent();
    assert(parent != NULL);
    newGenome = outAlignment->addLeafGenome(rootName, parent->getName(), 
                                            inAlignment->getBranchLength(
                                              parent->getName(), rootName));
  }
  assert(newGenome != NULL);
  
  vector<Sequence::Info> dimensions;
  getDimensions(genome, dimensions);
  newGenome->setDimensions(dimensions);

  // copy dna
  // copy top segments
  // copy bottom segments
  
  vector<string> childNames = inAlignment->getChildNames(rootName);
  for (size_t i = 0; i < childNames.size(); ++i)
  {
    extract(inAlignment, outAlignment, childNames[i]);
  }
}
