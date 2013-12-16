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

static void getDimensions(AlignmentConstPtr outAlignment, const Genome* genome,
                          vector<Sequence::Info>& dimensions);

static void copyGenome(const Genome* inGenome, Genome* outGenome);

static void extractTree(const AlignmentConstPtr inAlignment,
                        AlignmentPtr outAlignment, 
                        const string& rootName);

static void extract(const AlignmentConstPtr inAlignment,
                    AlignmentPtr outAlignment, const string& rootName);

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
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

    AlignmentPtr outAlignment = hdf5AlignmentInstance();
    outAlignment->setOptionsFromParser(optionsParser);
    outAlignment->createNew(outHalPath);

    if (outAlignment->getNumGenomes() != 0)
    {
      throw hal_exception("output hal alignmenet cannot be initialized");
    }
    
    if (rootName == "\"\"" || inAlignment->getNumGenomes() == 0)
    {
      rootName = inAlignment->getRootName();
    }

    extractTree(inAlignment, outAlignment, rootName);
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

void getDimensions(AlignmentConstPtr outAlignment, const Genome* genome, 
                   vector<Sequence::Info>& dimensions)
{
  assert(dimensions.size() == 0);
  SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
  SequenceIteratorConstPtr seqEndIt = genome->getSequenceEndIterator();

  bool root = outAlignment->getParentName(genome->getName()).empty();
  bool leaf = outAlignment->getChildNames(genome->getName()).empty();     
  
  for (; seqIt != seqEndIt; seqIt->toNext())
  {
    const Sequence* sequence = seqIt->getSequence();
    Sequence::Info info(sequence->getName(),
                        sequence->getSequenceLength(),
                        root ? 0 : sequence->getNumTopSegments(),
                        leaf ? 0 : sequence->getNumBottomSegments());
    dimensions.push_back(info);
  }
}

void copyGenome(const Genome* inGenome, Genome* outGenome)
{
  DNAIteratorConstPtr inDna = inGenome->getDNAIterator();
  DNAIteratorPtr outDna = outGenome->getDNAIterator();
  hal_size_t n = inGenome->getSequenceLength();
  assert(n == outGenome->getSequenceLength());
  for (; (hal_size_t)inDna->getArrayIndex() < n; inDna->toRight(), 
          outDna->toRight())
  {
    outDna->setChar(inDna->getChar());
  }

  TopSegmentIteratorConstPtr inTop = inGenome->getTopSegmentIterator();
  TopSegmentIteratorPtr outTop = outGenome->getTopSegmentIterator();
  n = outGenome->getNumTopSegments();
  assert(n == 0 || n == inGenome->getNumTopSegments());
  for (; (hal_size_t)inTop->getArrayIndex() < n; inTop->toRight(),
          outTop->toRight())
  {
    outTop->setCoordinates(inTop->getStartPosition(), inTop->getLength());
    outTop->setParentIndex(inTop->getParentIndex());
    outTop->setParentReversed(inTop->getParentReversed());
    outTop->setBottomParseIndex(inTop->getBottomParseIndex());
    outTop->setNextParalogyIndex(inTop->getNextParalogyIndex());
  }

  BottomSegmentIteratorConstPtr inBot = inGenome->getBottomSegmentIterator();
  BottomSegmentIteratorPtr outBot = outGenome->getBottomSegmentIterator();
  n = outGenome->getNumBottomSegments();
  assert(n == 0 || n == inGenome->getNumBottomSegments());
  hal_size_t nc = inGenome->getNumChildren();
  assert(nc == outGenome->getNumChildren());
  for (; (hal_size_t)inBot->getArrayIndex() < n; inBot->toRight(),
          outBot->toRight())
  {
    outBot->setCoordinates(inBot->getStartPosition(), inBot->getLength());
    for (hal_size_t child = 0; child < nc; ++child)
    {
      outBot->setChildIndex(child, inBot->getChildIndex(child));
      outBot->setChildReversed(child, inBot->getChildReversed(child));
    }
    if (outGenome->getAlignment()->getRootName() == outGenome->getName())
    {
      outBot->setTopParseIndex(NULL_INDEX);
    }
    else
    {
      outBot->setTopParseIndex(inBot->getTopParseIndex());
    }
  }

  const map<string, string>& meta = inGenome->getMetaData()->getMap();
  map<string, string>::const_iterator i = meta.begin();
  for (; i != meta.end(); ++i)
  {
    outGenome->getMetaData()->set(i->first, i->second);
  }
}

static void extractTree(const AlignmentConstPtr inAlignment,
                        AlignmentPtr outAlignment, 
                        const string& rootName)
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

  inAlignment->closeGenome(genome);
  outAlignment->closeGenome(newGenome);
  
  vector<string> childNames = inAlignment->getChildNames(rootName);
  for (size_t i = 0; i < childNames.size(); ++i)
  {
    extractTree(inAlignment, outAlignment, childNames[i]);
  }
}


void extract(const AlignmentConstPtr inAlignment,
             AlignmentPtr outAlignment, const string& rootName)
{
  const Genome* genome = inAlignment->openGenome(rootName);
  Genome* newGenome = outAlignment->openGenome(rootName);
  assert(newGenome != NULL);

  vector<Sequence::Info> dimensions;
  getDimensions(outAlignment, genome, dimensions);
  newGenome->setDimensions(dimensions);

  cout << "Extracting " << genome->getName() << endl;
  copyGenome(genome, newGenome);  

  inAlignment->closeGenome(genome);
  outAlignment->closeGenome(newGenome);
  
  vector<string> childNames = inAlignment->getChildNames(rootName);
  for (size_t i = 0; i < childNames.size(); ++i)
  {
    extract(inAlignment, outAlignment, childNames[i]);
  }
}
