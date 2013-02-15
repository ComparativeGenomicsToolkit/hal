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
static void printBranchPath(ostream& os, AlignmentConstPtr alignment, 
                            const vector<string>& genomeNames);
static void printBranches(ostream& os, AlignmentConstPtr alignment); 

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
  optionsParser->addOptionFlag("tree", "print only the NEWICK tree", false);
  optionsParser->addOptionFlag("branches", "print list of branches. "
                               "Each branch is specified by the child genome", 
                               false);
  optionsParser->addOption("path", "print branches on path between comma "
                           "separated list of genomes", "\"\"");
  string path;
  bool listGenomes;
  string sequencesFromGenome;
  string sequenceStatsFromGenome;
  string pathGenomes;
  bool tree;
  bool branches;
  try
  {
    optionsParser->parseOptions(argc, argv);
    path = optionsParser->getArgument<string>("halFile");
    listGenomes = optionsParser->getFlag("genomes");
    sequencesFromGenome = optionsParser->getOption<string>("sequences");
    sequenceStatsFromGenome = optionsParser->getOption<string>("sequenceStats");
    tree = optionsParser->getFlag("tree");
    pathGenomes = optionsParser->getOption<string>("path");
    branches = optionsParser->getFlag("branches");

    size_t optCount = listGenomes == true ? 1 : 0;
    if (sequencesFromGenome != "\"\"") ++optCount;
    if (tree == true) ++optCount;
    if (sequenceStatsFromGenome != "\"\"") ++optCount;
    if (pathGenomes != "\"\"") ++optCount;
    if (branches) ++ optCount;
    if (optCount > 1)
    {
      throw hal_exception("--genomes, --sequences, --tree, --path, --branches "
                          "and --sequenceStats " 
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
    else if (pathGenomes !=  "\"\"")
    {
      printBranchPath(cout, alignment, chopString(pathGenomes, ","));
    }
    else if (branches == true)
    {
      printBranches(cout, alignment);
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
  
  vector<const Genome*> outputVec;
  // if given two genomes, sort the output to be the actual path frmo the 
  // first to the second. 
  if (genomeNames.size() == 2)
  {
    outputVec.push_back(alignment->openGenome(genomeNames[0]));
    outputSet.erase(alignment->openGenome(genomeNames[0]));
    while (outputSet.empty() == false)
    {
      const Genome* cur = outputVec.back();
      set<const Genome*>::iterator i = outputSet.find(cur->getParent());
      if (i == outputSet.end())
      {
        for (size_t childIdx = 0; childIdx < cur->getNumChildren(); ++childIdx)
        {
          i = outputSet.find(cur->getChild(childIdx));
          if (i != outputSet.end())
          {
            break;
          }
        }
      }
      if (i != outputSet.end())
      {
        outputVec.push_back(*i);
        outputSet.erase(i);
      }
      else
      {
        throw hal_exception(string("error determining path from ") +
                            genomeNames[0] + " to " + genomeNames[1]);
      }
    }
  }
  else
  {
    outputVec.resize(outputSet.size());
    copy(outputSet.begin(), outputSet.end(), outputVec.begin());
  }
  
  for (vector<const Genome*>::const_iterator j = outputVec.begin(); 
       j != outputVec.end(); ++j)
  {
    const Genome* genome = *j;
    os << genome->getName() << " ";
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

