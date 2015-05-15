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
static void printGenomeMetaData(ostream &os, AlignmentConstPtr alignment,
                          const string &genomeName);
static void printChromSizes(ostream& os, AlignmentConstPtr alignment, 
                            const string& genomeName);
static void printPercentID(ostream& os, AlignmentConstPtr alignment,
                           const string& genomeName);
static void printCoverage(ostream& os, AlignmentConstPtr alignment,
                                 const string& genomeName);
static void printSegments(ostream& os, AlignmentConstPtr alignment,
                          const string& genomeName, bool top);
static void printAllCoverage(ostream& os, AlignmentConstPtr alignment);

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->setDescription("Retrieve basic statistics from a hal database");
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
  optionsParser->addOption("genomeMetaData", "print metadata for given genome, "
                           "one entry per line, tab-seperated.", "\"\"");
  optionsParser->addOption("chromSizes", "print the name and length of each"
                           " sequence in a given genome.  This is a subset"
                           " of the"
                           " information returned by --sequenceStats but is"
                           " useful because it is in the format used by"
                           " wigToBigWig", 
                           "\"\"");
  optionsParser->addOption("percentID",
                           "print % ID of a genome with all other genomes."
                           "Only non-duplicated and unambiguous sites are"
                           "considered",
                           "\"\"");
  optionsParser->addOption("coverage",
                           "print histogram of coverage of a genome with"
                           " all genomes", "\"\"");
  optionsParser->addOption("topSegments",
                           "print coordinates of all top segments of given"
                           " genome in BED format.", "\"\"");
  optionsParser->addOption("bottomSegments",
                           "print coordinates of all bottom segments of given"
                           " genome in BED format.", "\"\"");
  optionsParser->addOptionFlag("allCoverage",
                               "print histogram of coverage from all genomes to"
                               " all genomes", false);


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
  string genomeMetaData;
  string chromSizesFromGenome;
  string percentID;
  string coverage;
  string topSegments;
  string bottomSegments;
  bool allCoverage;
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
    genomeMetaData = optionsParser->getOption<string>("genomeMetaData");
    chromSizesFromGenome = optionsParser->getOption<string>("chromSizes");
    percentID = optionsParser->getOption<string>("percentID");
    coverage = optionsParser->getOption<string>("coverage");
    topSegments = optionsParser->getOption<string>("topSegments");
    bottomSegments = optionsParser->getOption<string>("bottomSegments");
    allCoverage = optionsParser->getFlag("allCoverage");

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
    if (genomeMetaData != "\"\"") ++optCount;
    if (chromSizesFromGenome != "\"\"") ++optCount;
    if (percentID != "\"\"") ++optCount;
    if (coverage != "\"\"") ++optCount;
    if (topSegments != "\"\"") ++optCount;
    if (bottomSegments != "\"\"") ++optCount;
    if (allCoverage) ++optCount;
    if (optCount > 1)
    {
      throw hal_exception("--genomes, --sequences, --tree, --span, --spanRoot, "
                          "--branches, --sequenceStats, --children, --parent, "
                          "--bedSequences, --root, --numSegments, --baseComp, "
                          "--genomeMetaData, --chromSizes, --percentID, "
                          "--coverage,  --topSegments, --bottomSegments, "
                          "--allCoverage "
                          "and --branchLength options are exclusive");
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
    else if (genomeMetaData != "\"\"")
    {
      printGenomeMetaData(cout, alignment, genomeMetaData);
    }
    else if (chromSizesFromGenome != "\"\"")
    {
      printChromSizes(cout, alignment, chromSizesFromGenome);
    }
    else if (percentID != "\"\"")
    {
      printPercentID(cout, alignment, percentID);
    }
    else if (coverage != "\"\"") {
      printCoverage(cout, alignment, coverage);
    }
    else if (topSegments != "\"\"") {
      printSegments(cout, alignment, topSegments, true);
    }
    else if (bottomSegments != "\"\"") {
      printSegments(cout, alignment, bottomSegments, false);
    } else if (allCoverage) {
      printAllCoverage(cout, alignment);
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

void printGenomeMetaData(ostream &os, AlignmentConstPtr alignment,
                         const string &genomeName)
{
  const Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
    throw hal_exception("Genome not found: " + genomeName);
  }
  const MetaData *metaData = genome->getMetaData();
  const map<string, string> metaDataMap = metaData->getMap();
  map<string, string>::const_iterator mapIt = metaDataMap.begin();
  for (; mapIt != metaDataMap.end(); mapIt++) {
    os << mapIt->first << '\t' << mapIt->second << '\n';
  }
  alignment->closeGenome(genome);
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

void printPercentID(ostream& os, AlignmentConstPtr alignment,
                    const string& genomeName)
{
  const Genome *refGenome = alignment->openGenome(genomeName);
  if (!refGenome) {
    throw hal_exception("Genome " + genomeName + " does not exist.");
  }

  ColumnIteratorPtr colIt = refGenome->getColumnIterator();
  // A bit sloppy, but a mapping from genome to (# identical bases, # aligned sites)
  // The # of aligned sites is necessary since a) not all sites are aligned and
  // b) we don't consider anything containing N's to be aligned.
  map<const Genome *, pair<hal_size_t *, hal_size_t *> > genomeStats;
  while(1) {
    // Get DNA for this site in reference
    DNAIteratorConstPtr refDnaIt = refGenome->getDNAIterator(colIt->getReferenceSequencePosition() + colIt->getReferenceSequence()->getStartPosition());
    char refDna = toupper(refDnaIt->getChar());
    
    const ColumnIterator::ColumnMap *cmap = colIt->getColumnMap();
    map <const Genome *, pair<hal_size_t *, hal_size_t *> > tempGenomeStats;
    for (ColumnIterator::ColumnMap::const_iterator colMapIt = cmap->begin();
         colMapIt != cmap->end(); colMapIt++) {
      if (colMapIt->second->empty()) {
        // There are empty entries in the column map.
        continue;
      }
      const Genome *genome = colMapIt->first->getGenome();
      const ColumnIterator::DNASet *dnaSet = colMapIt->second;
      assert(dnaSet->size() == 1);
      for (hal_size_t i = 0; i < dnaSet->size(); i++) {
        DNAIteratorConstPtr dnaIt = dnaSet->at(i);
        char otherDna = toupper(dnaIt->getChar());
        if (refDna != 'N' && otherDna != 'N') {
          if (!tempGenomeStats.count(genome)) {
            // initialize the map for this genome if necessary.
            tempGenomeStats[genome] = make_pair(new hal_size_t, new hal_size_t);
            *tempGenomeStats[genome].first = 0;
            *tempGenomeStats[genome].second = 0;
          }
          hal_size_t *tempNumID = tempGenomeStats[genome].first;
          hal_size_t *tempNumSites = tempGenomeStats[genome].second;
          if (refDna == otherDna) {
            (*tempNumID)++;
          }
          (*tempNumSites)++;
        }
      }
    }
    if (refDna != 'N' && *tempGenomeStats[refGenome].second == 1) {
      // If there isn't a duplication in the reference then we count
      // this column.
      for (map<const Genome *, pair<hal_size_t *, hal_size_t *> >::iterator it = tempGenomeStats.begin();
           it != tempGenomeStats.end(); it++) {
        const Genome *genome = it->first;
        if (!genomeStats.count(genome)) {
          // initialize the map for this genome if necessary.
          genomeStats[genome] = make_pair(new hal_size_t, new hal_size_t);
          *genomeStats[genome].first = 0;
          *genomeStats[genome].second = 0;
        }
        hal_size_t *tempNumID = it->second.first;
        hal_size_t *tempNumSites = it->second.second;
        if (*tempNumSites == 1) {
          // only count for the ID ratio if there's only 1 site for
          // this genome in the column.
          hal_size_t *numID = genomeStats[genome].first;
          hal_size_t *numSites = genomeStats[genome].second;
          assert(*tempNumID == 1 || *tempNumID == 0);
          assert(*tempNumSites == 1);
          if (*tempNumID == 1) {
            (*numID)++;
          }
          (*numSites)++;
        }
      }
    }
    // clean up temporary counts.
    for (map<const Genome *, pair<hal_size_t *, hal_size_t *> >::iterator it = tempGenomeStats.begin();
         it != tempGenomeStats.end(); it++) {
      delete it->second.first;
      delete it->second.second;
    }
    if (colIt->getReferenceSequencePosition() % 1000 == 0) {
      colIt->defragment();
    }
    if (colIt->lastColumn()) {
      // Break here--the column iterator will crash if we try to go further.
      break;
    }
    colIt->toRight();
  }
  os << "Genome, % ID, numID, numSites" << endl;
  for (map<const Genome *, pair<hal_size_t *, hal_size_t *> >::iterator statsIt = genomeStats.begin();
       statsIt != genomeStats.end(); statsIt++) {
    string name = statsIt->first->getName();
    hal_size_t numID = *statsIt->second.first;
    hal_size_t numSites = *statsIt->second.second;
    os << name << ", " << ((double) numID)/numSites << ", " << numID <<
      ", " << numSites << endl;
    delete statsIt->second.first;
    delete statsIt->second.second;
  }
}

void printCoverage(ostream& os, AlignmentConstPtr alignment,
                          const string& genomeName)
{
  const Genome *refGenome = alignment->openGenome(genomeName);
  if (!refGenome) {
    throw hal_exception("Genome " + genomeName + " does not exist.");
  }
  // Follow paralogies, but ignore ancestors.
  ColumnIteratorPtr colIt = refGenome->getColumnIterator(NULL, 0, 0,
                                                         NULL_INDEX, false,
                                                         true, false, true);
  map<const Genome *, vector<hal_size_t> *> histograms;
  while(1) {
    const ColumnIterator::ColumnMap *cmap = colIt->getColumnMap();
    // Temporary collecting of per-genome sites mapped, since it's
    // organized in the column map by sequence, not genome.
    map<const Genome *, hal_size_t> numSitesMapped;
    for (ColumnIterator::ColumnMap::const_iterator colMapIt = cmap->begin();
         colMapIt != cmap->end(); colMapIt++) {
      if (colMapIt->second->empty()) {
        // There are empty entries in the column map.
        continue;
      }
      const Genome *genome = colMapIt->first->getGenome();
      if (!numSitesMapped.count(genome)) {
        // Initialize map entry
        numSitesMapped[genome] = 0;
      }
      const ColumnIterator::DNASet *dnaSet = colMapIt->second;
      numSitesMapped[genome] = numSitesMapped[genome] + dnaSet->size();
    }
    for (map<const Genome *, hal_size_t>::const_iterator it = numSitesMapped.begin();
         it != numSitesMapped.end(); it++) {
      if (!histograms.count(it->first)) {
        // Initialize map
        histograms[it->first] = new vector<hal_size_t>;
      }
      vector<hal_size_t> *histogram = histograms[it->first];
      if (histogram->size() < it->second) {
        histogram->resize(it->second, 0);
      }
      for (hal_size_t i = 0; i < it->second; i++) {
        (*histogram)[i] = histogram->at(i) + numSitesMapped[refGenome];
      }
    }
    if (colIt->getReferenceSequencePosition() % 1000 == 0) {
      colIt->defragment();
    }
    if (colIt->lastColumn()) {
      // Break here--the column iterator will crash if we try to go further.
      break;
    }
    colIt->toRight();
  }
  hal_size_t maxHistLength = 0;
  for (map<const Genome *, vector<hal_size_t> *>::iterator histIt = histograms.begin();
       histIt != histograms.end(); histIt++) {
    vector <hal_size_t> *histogram = histIt->second;
    if (histogram->size() > maxHistLength) {
      maxHistLength = histogram->size();
    }
  }

  os << "Genome";
  for (hal_size_t i = 0; i < maxHistLength; i++) {
    os << ", sitesCovered" << i + 1 << "Times";
  }
  os << endl;
  for (map<const Genome *, vector<hal_size_t> *>::iterator histIt = histograms.begin();
       histIt != histograms.end(); histIt++) {
    string name = histIt->first->getName();
    os << name;
    vector <hal_size_t> *histogram = histIt->second;
    for(hal_size_t i = 0; i < maxHistLength; i++) {
      if (i < histogram->size()) {
        os << ", " << (double) histogram->at(i);
      } else {
        os << ", " << 0;
      }
    }
    os << endl;
  }
}

static void printSegments(ostream& os, AlignmentConstPtr alignment,
                          const string& genomeName, bool top)
{
  const Genome* genome = alignment->openGenome(genomeName);
  if (genome == NULL) 
  {
    throw hal_exception("Genome " + genomeName + " does not exist.");
  }
  hal_size_t numSegments = 0;
  SegmentIteratorConstPtr segment;
  if (top == true)
  {
    numSegments = genome->getNumTopSegments();
    if (numSegments > 0)
    {
      segment = genome->getTopSegmentIterator();
    }
  }
  else
  {
    numSegments = genome->getNumBottomSegments();
    if (numSegments > 0)
    {
      segment = genome->getBottomSegmentIterator();
    }
  }
  for (hal_size_t i = 0; i < numSegments; ++i)
  {
    const Sequence* sequence = segment->getSequence();
    os << sequence->getName() << '\t'
       << (segment->getStartPosition() - sequence->getStartPosition()) << '\t'
       << (segment->getEndPosition() + 1 - sequence->getStartPosition()) << '\n';
    segment->toRight();
  }
}

// Print coverage for all leaves vs. all leaves efficiently.
static void printAllCoverage(ostream& os, AlignmentConstPtr alignment)
{
  vector<string> leafNames = alignment->getLeafNamesBelow(alignment->getRootName());
  vector<const Genome *> leafGenomes;
  for (hal_size_t i = 0; i < leafNames.size(); i++) {
    const Genome *genome = alignment->openGenome(leafNames[i]);
    assert(genome != NULL);
    leafGenomes.push_back(genome);
  }
  ColumnIterator::VisitCache visitCache;
  map<pair<const Genome *, const Genome *>, vector<hal_size_t> *> histograms;
  for (hal_size_t i = 0; i < leafGenomes.size(); i++) {
    const Genome *genome = leafGenomes[i];
    // Follow paralogies, but ignore ancestors.
    ColumnIteratorPtr colIt = genome->getColumnIterator(NULL, 0, 0,
                                                        NULL_INDEX, false,
                                                        true, false, true);
    colIt->setVisitCache(&visitCache);
    // So that we don't accidentally visit the first column if it's
    // already been visited.
    colIt->toSite(0, genome->getSequenceLength() - 1);
    while(1) {
      const ColumnIterator::ColumnMap *cmap = colIt->getColumnMap();
      // Temporary collecting of per-genome sites mapped, since it's
      // organized in the column map by sequence, not genome.
      map<const Genome *, hal_size_t> numSitesMapped;
      for (ColumnIterator::ColumnMap::const_iterator colMapIt = cmap->begin();
           colMapIt != cmap->end(); colMapIt++) {
        if (colMapIt->second->empty()) {
          // There are empty entries in the column map.
          continue;
        }
        const Genome *genome = colMapIt->first->getGenome();
        if (!numSitesMapped.count(genome)) {
          // Initialize map entry
          numSitesMapped[genome] = 0;
        }
        const ColumnIterator::DNASet *dnaSet = colMapIt->second;
        numSitesMapped[genome] = numSitesMapped[genome] + dnaSet->size();
      }
      // O(n^2) in the number of genomes in the column -- doesn't seem
      // like there is a better way, since coverage isn't quite
      // symmetric.
      for (map<const Genome *, hal_size_t>::const_iterator it = numSitesMapped.begin();
           it != numSitesMapped.end(); it++) {
        for(map<const Genome *, hal_size_t>::const_iterator it2 = numSitesMapped.begin();
           it2 != numSitesMapped.end(); it2++) {
          pair<const Genome *, const Genome *> key = make_pair(it->first, it2->first);
          if (!histograms.count(key)) {
            // Initialize map
            histograms[key] = new vector<hal_size_t>;
          }
          vector<hal_size_t> *histogram = histograms[key];
          if (histogram->size() < it2->second) {
            histogram->resize(it2->second, 0);
          }
          for (hal_size_t i = 0; i < it2->second; i++) {
            (*histogram)[i] = histogram->at(i) + numSitesMapped[it->first];
          }
        }
      }
      if (colIt->getReferenceSequencePosition() % 1000 == 0) {
        colIt->defragment();
      }
      if (colIt->lastColumn()) {
        // Break here--the column iterator will crash if we try to go further.
        break;
      }
      colIt->toRight();
    }
    // Copy over the updated visit cache information so we can supply it to the next genome.
    visitCache.clear();
    ColumnIterator::VisitCache *newVisitCache = colIt->getVisitCache();
    for(ColumnIterator::VisitCache::iterator it = newVisitCache->begin();
        it != newVisitCache->end(); it++) {
      visitCache[it->first] = new PositionCache(*it->second);
    }
  }

  hal_size_t maxHistLength = 0;
  for (map<pair<const Genome *, const Genome *>, vector<hal_size_t> *>::iterator histIt = histograms.begin();
       histIt != histograms.end(); histIt++) {
    vector <hal_size_t> *histogram = histIt->second;
    if (histogram->size() > maxHistLength) {
      maxHistLength = histogram->size();
    }
  }

  os << "FromGenome, ToGenome";
  for (hal_size_t i = 0; i < maxHistLength; i++) {
    os << ", sitesCovered" << i + 1 << "Times";
  }
  os << endl;
  for (map<pair<const Genome *, const Genome *>, vector<hal_size_t> *>::iterator histIt = histograms.begin();
       histIt != histograms.end(); histIt++) {
    string fromName = histIt->first.second->getName();
    string toName = histIt->first.first->getName();
    os << fromName;
    os << ", " << toName;
    vector <hal_size_t> *histogram = histIt->second;
    for(hal_size_t i = 0; i < maxHistLength; i++) {
      if (i < histogram->size()) {
        os << ", " << (double) histogram->at(i);
      } else {
        os << ", " << 0;
      }
    }
    os << endl;
  }
}
