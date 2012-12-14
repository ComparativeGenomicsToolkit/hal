/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include "halStats.h"

using namespace std;
using namespace hal;

static void printGenomes(ostream& os, AlignmentConstPtr alignment);
static void printSequences(ostream& os, AlignmentConstPtr alignment, 
                          const string& genomeName);
static void printSequenceStats(ostream& os, AlignmentConstPtr alignment, 
                               const string& genomeName);
static void printBedSequenceStats(ostream& os, AlignmentConstPtr alignment, 
                                  const string& genomeName);
static void printBranchPath(ostream& os, AlignmentConstPtr alignment, 
                            const vector<string>& genomeNames);
static void printBranches(ostream& os, AlignmentConstPtr alignment);
static void printChildren(ostream& os, AlignmentConstPtr alignment, 
                          const string& genomeName);
static void printParent(ostream& os, AlignmentConstPtr alignment, 
                        const string& genomeName);
static void printRootName(ostream& os, AlignmentConstPtr alignment);


int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->setDescription("Rertrieve basic statics from a hal database");
  optionsParser->addArgument("halFile", "path to hal file to analyze");
  optionsParser->addOptionFlag("genomes", "print only a list of genomes "
                               "in alignment", false);
  optionsParser->addOption("sequences", "print list of sequences in given "
                           "genome", "\"\"");
  optionsParser->addOption("sequenceStats", "print stats for each sequence in "
                           "given genome", "\"\"");
  optionsParser->addOption("bedSequences", "print sequences of given genome "
                           "in bed format",
                           "\"\"");
  optionsParser->addOptionFlag("tree", "print only the NEWICK tree", false);
  optionsParser->addOptionFlag("branches", "print list of branches. "
                               "Each branch is specified by the child genome", 
                               false);
  optionsParser->addOption("path", "print branches on path between comma "
                           "separated list of genomes", "\"\"");
  optionsParser->addOption("children", "print names of children of given "
                           "genome", "\"\"");
  optionsParser->addOptionFlag("root", "print root genome name", false);
  optionsParser->addOption("parent", "print name of parent of given genome",
                           "\"\"");

  string path;
  bool listGenomes;
  string sequencesFromGenome;
  string sequenceStatsFromGenome;
  string bedSequencesFromGenome;
  string pathGenomes;
  bool tree;
  bool branches;
  string childrenFromGenome;
  string parentFromGenome;
  bool printRoot;
  try
  {
    optionsParser->parseOptions(argc, argv);
    path = optionsParser->getArgument<string>("halFile");
    listGenomes = optionsParser->getFlag("genomes");
    sequencesFromGenome = optionsParser->getOption<string>("sequences");
    sequenceStatsFromGenome = optionsParser->getOption<string>("sequenceStats");
    bedSequencesFromGenome = optionsParser->getOption<string>("bedSequences");
    tree = optionsParser->getFlag("tree");
    pathGenomes = optionsParser->getOption<string>("path");
    branches = optionsParser->getFlag("branches");
    childrenFromGenome = optionsParser->getOption<string>("children");
    parentFromGenome = optionsParser->getOption<string>("parent");
    printRoot = optionsParser->getFlag("root");

    size_t optCount = listGenomes == true ? 1 : 0;
    if (sequencesFromGenome != "\"\"") ++optCount;
    if (tree == true) ++optCount;
    if (sequenceStatsFromGenome != "\"\"") ++optCount;
    if (bedSequencesFromGenome != "\"\"") ++optCount;
    if (pathGenomes != "\"\"") ++optCount;
    if (branches) ++ optCount;
    if (childrenFromGenome != "\"\"") ++optCount;
    if (parentFromGenome != "\"\"") ++optCount;
    if (printRoot) ++optCount;
    if (optCount > 1)
    {
      throw hal_exception("--genomes, --sequences, --tree, --path, --branches, "
                          "--sequenceStats, --children, --parent, "
                          "--bedSequences abd --root " 
                          "options are mutually exclusive");
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
    AlignmentConstPtr alignment = openHalAlignmentReadOnly(path, optionsParser);

    if (listGenomes == true && alignment->getNumGenomes() > 0)
    {
      printGenomes(cout, alignment);
    }
    else if (sequencesFromGenome != "\"\"")
    {
      printSequences(cout, alignment, sequencesFromGenome);
    }
    else if (tree == true)
    {
      cout << alignment->getNewickTree() << endl;
    }
    else if (sequenceStatsFromGenome != "\"\"")
    {
      printSequenceStats(cout, alignment, sequenceStatsFromGenome);
    }
    else if (bedSequencesFromGenome != "\"\"")
    {
      printBedSequenceStats(cout, alignment, bedSequencesFromGenome);
    }
    else if (pathGenomes !=  "\"\"")
    {
      printBranchPath(cout, alignment, chopString(pathGenomes, ","));
    }
    else if (branches == true)
    {
      printBranches(cout, alignment);
    }
    else if (childrenFromGenome != "\"\"")
    {
      printChildren(cout, alignment, childrenFromGenome);
    }
    else if (parentFromGenome != "\"\"")
    {
      printParent(cout, alignment, parentFromGenome);
    }
    else if (printRoot == true)
    {
      printRootName(cout, alignment);
    }
    else
    {
      HalStats halStats(alignment);
      cout << endl << "hal v" << alignment->getVersion() << "\n" << halStats;
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

void printGenomes(ostream& os, AlignmentConstPtr alignment)
{
  const Genome* root = alignment->openGenome(alignment->getRootName());
  set<const Genome*> genomes;
  getGenomesInSubTree(root, genomes);
  genomes.insert(root);
  for (set<const Genome*>::iterator i = genomes.begin(); i != genomes.end();
       ++i)
  {
    if (i != genomes.begin())
    {
      os << ",";
    }
    os << (*i)->getName();
  }
  os << endl;      
}

void printSequences(ostream& os, AlignmentConstPtr alignment, 
                   const string& genomeName)
{
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome ") + genomeName + " not found.");
  }
  if (genome->getNumSequences() > 0)
  {
    SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
    SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
    for (; !seqIt->equals(seqEnd); seqIt->toNext())
    {
      if (!seqIt->equals(genome->getSequenceIterator()))
      {
        os << ",";
      }
      os << seqIt->getSequence()->getName();
    }
  }
  os << endl;
}

void printSequenceStats(ostream& os, AlignmentConstPtr alignment, 
                        const string& genomeName)
{
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome ") + genomeName + " not found.");
  }
  if (genome->getNumSequences() > 0)
  {
    SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
    SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();
    os << "SequenceName, Length, NumTopSegments, NumBottomSegments" << endl;

    for (; !seqIt->equals(seqEnd); seqIt->toNext())
    {
      os << seqIt->getSequence()->getName() << ", "
         << seqIt->getSequence()->getSequenceLength() << ", "
         << seqIt->getSequence()->getNumTopSegments() << ", "
         << seqIt->getSequence()->getNumBottomSegments() << "\n";
    }
  }
  os << endl;
}

void printBedSequenceStats(ostream& os, AlignmentConstPtr alignment, 
                           const string& genomeName)
{
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome ") + genomeName + " not found.");
  }
  if (genome->getNumSequences() > 0)
  {
    SequenceIteratorConstPtr seqIt = genome->getSequenceIterator();
    SequenceIteratorConstPtr seqEnd = genome->getSequenceEndIterator();

    for (; !seqIt->equals(seqEnd); seqIt->toNext())
    {
      os << seqIt->getSequence()->getName() << "\t"
         << 0 << "\t"
         << seqIt->getSequence()->getSequenceLength() << "\n";
    }
  }
  os << endl;
}

static void printBranchPath(ostream& os, AlignmentConstPtr alignment, 
                            const vector<string>& genomeNames)
{
  set<const Genome*> inputSet;
  for (size_t i = 0; i < genomeNames.size(); ++i)
  {
    const Genome* genome = alignment->openGenome(genomeNames[i]);
    if (genome == NULL)
    {
      throw hal_exception(string("Genome ") + genomeNames[i] + " not found");
    }
    inputSet.insert(genome);
  }
  set<const Genome*> outputSet;
  getGenomesInSpanningTree(inputSet, outputSet);
  
  for (set<const Genome*>::const_iterator j = outputSet.begin(); 
       j != outputSet.end(); ++j)
  {
    const Genome* genome = *j;
    if (genome->getParent() != NULL && 
        outputSet.find(genome->getParent()) != outputSet.end())
    {
      os << genome->getName() << " ";
    }
  }
  os << endl;
}

static void printBranches(ostream& os, AlignmentConstPtr alignment)
{
  const Genome* root = alignment->openGenome(alignment->getRootName());
  set<const Genome*> genomes;
  getGenomesInSubTree(root, genomes);
  genomes.insert(root);
  bool first = true;
  for (set<const Genome*>::iterator i = genomes.begin(); i != genomes.end();
       ++i)
  {
    if ((*i)->getParent() != NULL)
    {
      if (!first)
      {
        os << " ";
      }
      else
      {
        first = false;
      }
      os << (*i)->getName();
    }
  }
  os << endl;      
}

void printChildren(ostream& os, AlignmentConstPtr alignment, 
                   const string& genomeName)
{
  vector<string> children = alignment->getChildNames(genomeName);
  for (size_t i = 0; i < children.size(); ++i)
  {
    os << children[i];
    if (i != children.size() - 1)
    {
      os << " ";
    }
  }
  os << endl;
}

void printParent(ostream& os, AlignmentConstPtr alignment, 
                        const string& genomeName)
{
  if (genomeName != alignment->getRootName())
  {
    os << alignment->getParentName(genomeName) << endl;
  }
}

void printRootName(ostream& os, AlignmentConstPtr alignment)
{
  os << alignment->getRootName() << endl;
}
