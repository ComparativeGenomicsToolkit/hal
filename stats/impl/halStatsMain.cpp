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
                            const vector<string>& genomeNames, bool keepRoot);
static void printBranches(ostream& os, AlignmentConstPtr alignment);
static void printChildren(ostream& os, AlignmentConstPtr alignment, 
                          const string& genomeName);
static void printParent(ostream& os, AlignmentConstPtr alignment, 
                        const string& genomeName);
static void printRootName(ostream& os, AlignmentConstPtr alignment);
static void printBranchLength(ostream& os, AlignmentConstPtr alignment, 
                              const string& genomeName);
static void printBranches(ostream& os, AlignmentConstPtr alignment); 
static void printNumSegments(ostream& os, AlignmentConstPtr alignment,
                             const string& genomeName); 
static void printBaseComp(ostream& os, AlignmentConstPtr alignment, 
                          const string& baseCompPair);
static void printChromSizes(ostream& os, AlignmentConstPtr alignment, 
                            const string& genomeName);


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
  optionsParser->addOption("span", "print branches on path (or spanning tree) "
                           "between comma "
                           "separated list of genomes", "\"\"");
  optionsParser->addOption("spanRoot", "print genomes on path" 
                           "(or spanning tree) between comma "
                           "separated list of genomes.  Different from --span"
                           "only in that the spanning tree root is also "
                           "given", "\"\"");
  optionsParser->addOption("children", "print names of children of given "
                           "genome", "\"\"");
  optionsParser->addOptionFlag("root", "print root genome name", false);
  optionsParser->addOption("parent", "print name of parent of given genome",
                           "\"\"");
  optionsParser->addOption("branchLength", "print branch length between "
                           "given genome and its parent in the tree",
                           "\"\"");
  optionsParser->addOption("numSegments", "print numTopSegments "
                           "numBottomSegments for given genome.",
                           "\"\"");
  optionsParser->addOption("baseComp", "print base composition for given "
                           "genome by sampling every step bases. Parameter "
                           "value is of the form genome,step.  Ex: "
                           "--baseComp human,1000.  The ouptut is of the form "
                           "fraction_of_As fraction_of_Gs fraction_of_Cs "
                           "fraction_of_Ts.", 
                           "\"\"");
  optionsParser->addOption("chromSizes", "print the name and length of each"
                           " sequence in a given genome.  This is a subset"
                           " of the"
                           " information returned by --sequenceStats but is"
                           " useful because it is in the format used by"
                           " wigToBigWig", 
                           "\"\"");


  string path;
  bool listGenomes;
  string sequencesFromGenome;
  string sequenceStatsFromGenome;
  string bedSequencesFromGenome;
  string spanGenomes;
  string spanRootGenomes;
  bool tree;
  bool branches;
  string childrenFromGenome;
  string parentFromGenome;
  bool printRoot;
  string nameForBL;
  string numSegmentsGenome;
  string baseCompPair;
  string chromSizesFromGenome;
  try
  {
    optionsParser->parseOptions(argc, argv);
    path = optionsParser->getArgument<string>("halFile");
    listGenomes = optionsParser->getFlag("genomes");
    sequencesFromGenome = optionsParser->getOption<string>("sequences");
    sequenceStatsFromGenome = optionsParser->getOption<string>("sequenceStats");
    bedSequencesFromGenome = optionsParser->getOption<string>("bedSequences");
    tree = optionsParser->getFlag("tree");
    spanGenomes = optionsParser->getOption<string>("span");
    spanRootGenomes = optionsParser->getOption<string>("spanRoot");
    branches = optionsParser->getFlag("branches");
    childrenFromGenome = optionsParser->getOption<string>("children");
    parentFromGenome = optionsParser->getOption<string>("parent");
    printRoot = optionsParser->getFlag("root");
    nameForBL = optionsParser->getOption<string>("branchLength");
    numSegmentsGenome = optionsParser->getOption<string>("numSegments");
    baseCompPair = optionsParser->getOption<string>("baseComp");
    chromSizesFromGenome = optionsParser->getOption<string>("chromSizes");

    size_t optCount = listGenomes == true ? 1 : 0;
    if (sequencesFromGenome != "\"\"") ++optCount;
    if (tree == true) ++optCount;
    if (sequenceStatsFromGenome != "\"\"") ++optCount;
    if (bedSequencesFromGenome != "\"\"") ++optCount;
    if (spanGenomes != "\"\"") ++optCount;
    if (spanRootGenomes != "\"\"") ++optCount;
    if (branches) ++ optCount;
    if (childrenFromGenome != "\"\"") ++optCount;
    if (parentFromGenome != "\"\"") ++optCount;
    if (printRoot) ++optCount;
    if (nameForBL != "\"\"") ++optCount;
    if (numSegmentsGenome != "\"\"") ++optCount;
    if (baseCompPair != "\"\"") ++optCount;
    if (chromSizesFromGenome != "\"\"") ++optCount;
    if (optCount > 1)
    {
      throw hal_exception("--genomes, --sequences, --tree, --span, --spanRoot, "
                          "--branches, --sequenceStats, --children, --parent, "
                          "--bedSequences, --root, --numSegments, --baseComp, "
                          "--chromSizes "
                          "and --branchLength options are exclusive" );
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
    else if (spanGenomes !=  "\"\"")
    {
      printBranchPath(cout, alignment, chopString(spanGenomes, ","), false);
    }
    else if (spanRootGenomes !=  "\"\"")
    {
      printBranchPath(cout, alignment, chopString(spanRootGenomes, ","), true);
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
    else if (nameForBL != "\"\"")
    {
      printBranchLength(cout, alignment, nameForBL);
    }
    else if (numSegmentsGenome != "\"\"")
    {
      printNumSegments(cout, alignment, numSegmentsGenome);
    }
    else if (baseCompPair != "\"\"")
    {
      printBaseComp(cout, alignment, baseCompPair);
    }
    else if (chromSizesFromGenome != "\"\"")
    {
      printChromSizes(cout, alignment, chromSizesFromGenome);
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
    set<const Genome*>::iterator next = i;
    ++next;   
    os << (*i)->getName();
    if (next != genomes.end())
    {
      os << " ";
    }
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
                            const vector<string>& genomeNames, 
                            bool keepRoot)
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
    set<const Genome*> visitSet(outputSet);
    outputVec.push_back(alignment->openGenome(genomeNames[0]));
    visitSet.erase(alignment->openGenome(genomeNames[0]));
    while (outputVec.back()->getName() != genomeNames[1])
    {
      const Genome* cur = outputVec.back();
      set<const Genome*>::iterator i = visitSet.find(cur->getParent());
      if (i == visitSet.end())
      {
        for (size_t childIdx = 0; childIdx < cur->getNumChildren(); ++childIdx)
        {
          i = visitSet.find(cur->getChild(childIdx));
          if (i != visitSet.end())
          {
            break;
          }
        }
      }
      if (i != visitSet.end())
      {
        outputVec.push_back(*i);
        visitSet.erase(i);
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
    if (keepRoot == true || 
        (genome->getParent() != NULL &&  
         outputSet.find(genome->getParent()) != outputSet.end()))
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

void printBranchLength(ostream& os, AlignmentConstPtr alignment, 
                       const string& genomeName)
{
  if (genomeName != alignment->getRootName())
  {
    string parentName = alignment->getParentName(genomeName);
    os << alignment->getBranchLength(parentName, genomeName) << endl;
  }
}

void printNumSegments(ostream& os, AlignmentConstPtr alignment, 
                      const string& genomeName)
{
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome ") + genomeName + " not found.");
  }
  os << genome->getNumTopSegments() << " " << genome->getNumBottomSegments()
     << endl;
}

void printBaseComp(ostream& os, AlignmentConstPtr alignment, 
                   const string& baseCompPair)
{
  string genomeName;
  hal_size_t step = 0;
  vector<string> tokens = chopString(baseCompPair, ",");
  if (tokens.size() == 2)
  {
    genomeName = tokens[0];
    stringstream ss(tokens[1]);
    ss >> step;
  }
  if (step == 0)
  {
    stringstream ss;
    ss << "Invalid value for --baseComp: " << baseCompPair << ".  Must be of"
       << " format genomeName,step";
    throw hal_exception(ss.str());
  }
      
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL)
  {
    throw hal_exception(string("Genome ") + genomeName + " not found.");
  }
  hal_size_t numA = 0;
  hal_size_t numC = 0;
  hal_size_t numG = 0;
  hal_size_t numT = 0;

  hal_size_t len = genome->getSequenceLength();
  if (step >= len)
  {
    step = len - 1;
  }

  DNAIteratorConstPtr dna = genome->getDNAIterator();
  for (hal_size_t i = 0; i < len; i += step)
  {
    dna->jumpTo(i);
    switch (dna->getChar())
    {
    case 'a':
    case 'A':
      ++numA;
      break;
    case 'c':
    case 'C':
      ++numC;
      break;
    case 'g':
    case 'G':
      ++numG;
      break;
    case 't':
    case 'T':
      ++numT;
      break;
    default:
      break;
    }
  } 
  
  double total = numA + numC + numG + numT;
  os << (double)numA / total << '\t'
     << (double)numC / total << '\t'
     << (double)numG / total << '\t'
     << (double)numT / total << '\n';
}

void printChromSizes(ostream& os, AlignmentConstPtr alignment, 
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
      os << seqIt->getSequence()->getName() << '\t'
         << seqIt->getSequence()->getSequenceLength() << '\n';
    }
  }
}
