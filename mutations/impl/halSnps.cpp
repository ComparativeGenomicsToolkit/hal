/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "hal.h"

using namespace std;
using namespace hal;

static void countSnps(const Genome* refGenome,
                      set<const Genome *> &targetGenomes,
                      hal_index_t start, hal_size_t length, bool doDupes,
                      ostream& refTsvStream,
                      map<const Genome *, hal_size_t> *numSnps,
                      map<const Genome *, hal_size_t> *numOrthologousPairs,
                      bool unique, hal_size_t minSpeciesForSnp);

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("refGenome",
                             "name of reference genome.");
  optionsParser->addArgument("targetGenomes",
                             "names of query genomes, comma-separated.");
  optionsParser->addOption("tsv",
                           "path of tsv file to write snps to in reference "
                           "genome coordinates, and containing the base "
                           "assignments for each genome",
                           "\"\"");
  optionsParser->addOption("noDupes",
                           "do not consider paralogies while mapping",
                           "\"\"");
  optionsParser->addOption("refSequence",
                           "name of reference sequence within reference genome"
                           " (all sequences if empty)",
                           "\"\"");
  optionsParser->addOption("start",
                           "coordinate within reference genome (or sequence"
                           " if specified) to start at",
                           0);
  optionsParser->addOption("length",
                           "length of the reference genome (or sequence"
                           " if specified) to convert.  If set to 0,"
                           " the entire thing is converted",
                           0);
  optionsParser->addOption("minSpeciesForSnp", "minimum number of "
                           "species that must have a different base than the "
                           "reference for a SNP to be reported in the TSV", 1);
  optionsParser->addOptionFlag("unique",
                               "Whether to ignore columns that are not "
                               "canonical on the reference genome", false);
  optionsParser->setDescription("Count snps between orthologous positions "
                                "in multiple genomes.  Outputs "
                                "targetGenome totalSnps totalCleanOrthologousPairs");
  return optionsParser;
}

int main(int argc, char** argv)
{
  CLParserPtr optionsParser = initParser();

  string halPath;
  string refGenomeName;
  string targetGenomesString;
  string tsvPath;
  bool noDupes;
  string refSequenceName;
  hal_index_t start;
  hal_size_t length;
  bool unique;
  hal_size_t minSpeciesForSnp;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    targetGenomesString = optionsParser->getArgument<string>("targetGenomes");
    tsvPath = optionsParser->getOption<string>("tsv");
    noDupes = optionsParser->getOption<bool>("noDupes");
    refSequenceName = optionsParser->getOption<string>("refSequence");
    start = optionsParser->getOption<hal_index_t>("start");
    length = optionsParser->getOption<hal_size_t>("length");
    minSpeciesForSnp = optionsParser->getOption<hal_size_t>("minSpeciesForSnp");
    unique = optionsParser->getFlag("unique");
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
      throw hal_exception("hal alignment is empty");
    }

    const Genome* refGenome = NULL;
    refGenome = alignment->openGenome(refGenomeName);
    if (refGenome == NULL)
    {
      throw hal_exception(string("Reference genome, ") + refGenomeName +
                          ", not found in alignment");
    }

    vector<string> targetGenomeNames = chopString(targetGenomesString, ",");
    set<const Genome *> targetGenomes;
    for (hal_size_t i = 0; i < targetGenomeNames.size(); i++)
    {
      const Genome *genome = alignment->openGenome(targetGenomeNames[i]);
      if (genome == NULL)
      {
        throw hal_exception("Target genome " + targetGenomeNames[i]
                            + " not found in alignment.");
      }
      targetGenomes.insert(genome);
    }

    if (start + length >= refGenome->getSequenceLength())
    {
      throw hal_exception(string("Invalid range for ") + refGenomeName);
    }

    const Sequence* refSequence = NULL;
    if (refSequenceName != "\"\"")
    {
      refSequence = refGenome->getSequence(refSequenceName);
      if (refSequence == NULL)
      {
        throw hal_exception(string("Reference sequence, ") + refSequenceName +
                            ", not found in reference genome, " +
                            refGenome->getName());
      }
      if (start + length >= refSequence->getSequenceLength())
      {
        throw hal_exception(string("Invalid range for ") + refSequenceName);
      }
      if (length == 0)
      {
        length = refSequence->getSequenceLength() - start;
      }
      start += refSequence->getStartPosition();
    }

    if (length == 0)
    {
      assert(refSequence == NULL);
      length = refGenome->getSequenceLength() - start;
    }

    ofstream refTsvStream;
    if (tsvPath != "\"\"")
    {
      refTsvStream.open(tsvPath.c_str());
      if (!refTsvStream)
      {
        throw hal_exception("Error opening " + tsvPath);
      }
    }
    // build and initialize the snps/orthologous pairs maps per
    // genome.
    map<const Genome *, hal_size_t> numSnps;
    map<const Genome *, hal_size_t> numOrthologousPairs;

    for (set<const Genome *>::const_iterator i = targetGenomes.begin();
         i != targetGenomes.end(); i++)
    {
      numSnps[*i] = 0;
      numOrthologousPairs[*i] = 0;
    }

    countSnps(refGenome, targetGenomes, start, length, !noDupes,
              refTsvStream, &numSnps, &numOrthologousPairs,
              unique, minSpeciesForSnp);

    assert(numSnps.size() == numOrthologousPairs.size());

    for (map<const Genome *, hal_size_t>::const_iterator snpIt = numSnps.begin();
         snpIt !=  numSnps.end(); snpIt++)
    {
      const Genome *genome = snpIt->first;
      hal_size_t numSnps = snpIt->second;
      hal_size_t orthologousPairs = numOrthologousPairs[genome];
      cout << genome->getName() << " " << numSnps << " "
           << orthologousPairs << endl;
      
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

// Recursively find nodes that are from the reference genome in the
// tree and add them to the refNodes set.
static void getReferenceNodes_R(stTree *colTree, const Genome *refGenome,
                                set<stTree *> *refNodes)
{
  for (int64_t i = 0; i < stTree_getChildNumber(colTree); i++)
  {
    getReferenceNodes_R(stTree_getChild(colTree, i), refGenome, refNodes);
  }

  DNAIteratorConstPtr *dnaIt = (DNAIteratorConstPtr *) stTree_getClientData(colTree);
  if ((*dnaIt)->getGenome() == refGenome) {
    assert((refGenome->getNumChildren() != 0) ^
           (stTree_getChildNumber(colTree) == 0));
    refNodes->insert(colTree);
  }
}

// copy of stTree_getMRCA since phylogeny branch isn't merged in yet.
static stTree *stTree_getMRCA_tmp(stTree *node1, stTree *node2) {
    // Find all of node 1's parents (inclusive of node 1)
    stSet *parents = stSet_construct();
    stTree *curNode = node1;
    do {
        stSet_insert(parents, curNode);
    } while ((curNode = stTree_getParent(curNode)) != NULL);

    // Find the first parent of node 2 that is a parent of node 1
    stTree *ret = NULL;
    curNode = node2;
    do {
        if (stSet_search(parents, curNode) != NULL) {
            ret = curNode;
            break;
        }
    } while ((curNode = stTree_getParent(curNode)) != NULL);

    stSet_destruct(parents);
    return ret;
}

static void addSubtreeToOrthologSet(stTree *tree, set <const Genome *> &targetGenomes, set<DNAIteratorConstPtr *> *orthologSet)
{
  DNAIteratorConstPtr *dnaIt = (DNAIteratorConstPtr *) stTree_getClientData(tree);
  if (targetGenomes.count((*dnaIt)->getGenome()))
  {
    orthologSet->insert(dnaIt);
  }
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++)
  {
    addSubtreeToOrthologSet(stTree_getChild(tree, i), targetGenomes, orthologSet);
  }
}

// Remove any ptrs to DNAIteratorConstPtrs that belong to the same
// genome as another entry in the set.
static void removeDuplicatedGenomes(set<DNAIteratorConstPtr *> *orthologSet)
{
  set<const Genome *> seen, genomesToRemove;

  // Find genomes that have more than one DNAIterator in the set.
  for (set<DNAIteratorConstPtr *>::const_iterator i = orthologSet->begin();
       i != orthologSet->end(); i++) {
    const Genome *genome = (**i)->getGenome();
    if (seen.count(genome))
    {
      genomesToRemove.insert(genome);
    }
    seen.insert(genome);
  }

  // Go through the set again and remove any genomes that have more
  // than one element.
  for (set<DNAIteratorConstPtr *>::const_iterator i = orthologSet->begin();
       i != orthologSet->end(); ) {
    const Genome *genome = (**i)->getGenome();
    if (genomesToRemove.count(genome)) {
      // remove the current element, but advance the iterator first so
      // it doesn't become invalidated by erase()
      orthologSet->erase(i++);
    } else {
      i++;
    }
  }
}

// Get a mapping from ref genome base -> orthologous non-reference
// bases from a column tree. Only clear orthologs are added.
static void getOrthologs(stTree *colTree, const Genome *refGenome,
                         set<const Genome *> &targetGenomes,
                         map<DNAIteratorConstPtr *, set<DNAIteratorConstPtr *> *> *orthologs)
{
//  printf("%s\n", stTree_getNewickTreeString(colTree));
  set<stTree *> refNodes;
  getReferenceNodes_R(colTree, refGenome, &refNodes);

  // now find all orthologous bases. Additionally we get rid of any
  // duplications that have happened after diverging from the MRCA.

  // Get the set of coalescences of all pairs of the ref nodes.
  set<stTree *>refCoalescences;
  for (set<stTree *>::const_iterator refNodeIt = refNodes.begin(); refNodeIt != refNodes.end(); refNodeIt++)
  {
    for (set<stTree *>::const_iterator refNodeIt2 = refNodeIt; refNodeIt2 != refNodes.end(); refNodeIt2++) {
      refCoalescences.insert(stTree_getMRCA_tmp(*refNodeIt, *refNodeIt2));
    }

    // Use those coalescences as "stops" and traverse up the tree from
    // the ref nodes. The maximal subtrees below a ref node containing
    // one ref node will contain the orthologs.
    DNAIteratorConstPtr *refDnaIt = (DNAIteratorConstPtr *) stTree_getClientData(*refNodeIt);
    set<DNAIteratorConstPtr *> *orthologSet = new set<DNAIteratorConstPtr *>();
    (*orthologs)[refDnaIt] = orthologSet;
    stTree *curNode = *refNodeIt;
    while (stTree_getParent(curNode) != NULL && refCoalescences.count(stTree_getParent(curNode)) == 0)
    {
      curNode = stTree_getParent(curNode);
    }

    // OK, we've found the root of the subtree containing the
    // orthologs, add them to the set
    addSubtreeToOrthologSet(curNode, targetGenomes, orthologSet);

    // Get rid of cases where there is more than one ortholog per
    // genome, i.e. there has been a duplication since the MRCA of the
    // reference node and a target node.
    removeDuplicatedGenomes(orthologSet);
  }
}

static void countSnps(const Genome* refGenome,
                      set<const Genome *> &targetGenomes,
                      hal_index_t start, hal_size_t length, bool doDupes,
                      ostream& refTsvStream,
                      map<const Genome *, hal_size_t> *numSnps,
                      map<const Genome *, hal_size_t> *numOrthologousPairs,
                      bool unique, hal_size_t minSpeciesForSnp)
{

  // A bit annoying, but we have to build a map from genome to tsv
  // field position so that the ordering is consistent in the tsv.
  // We also print out the tsv header while we're at it.
  map<const Genome *, hal_size_t> fieldForGenome;
  if (refTsvStream) {
    refTsvStream << "refSequence\trefPosition\t";
    fieldForGenome[refGenome] = 2;
    refTsvStream << refGenome->getName();
    int64_t fieldNum = 3;
    for (set<const Genome *>::const_iterator i = targetGenomes.begin();
         i != targetGenomes.end(); i++) {
      refTsvStream << "\t" << (*i)->getName();
      fieldForGenome[*i] = fieldNum;
      fieldNum++;
    }
    refTsvStream << endl;
  }

  hal_index_t lastPos = start + length - 1;

  ColumnIteratorConstPtr colIt = refGenome->getColumnIterator(&targetGenomes,
                                                              0,
                                                              start,
                                                              lastPos,
                                                              !doDupes,
                                                              true);
  while (1)
  {
    if (unique && !colIt->isCanonicalOnRef()) {
      // This column isn't unique (if we iterate over the reference
      // segments separately, we will have visited this column
      // already).
      continue;
    }
    stTree *colTree = colIt->getTree();
    map<DNAIteratorConstPtr *, set<DNAIteratorConstPtr *> *> orthologs;
    getOrthologs(colTree, refGenome, targetGenomes, &orthologs);

    // Now that we have the set of reference bases and their
    // orthologs, just call SNPs.
    for (map<DNAIteratorConstPtr *, set<DNAIteratorConstPtr *> *>::const_iterator orthologsIt = orthologs.begin(); orthologsIt != orthologs.end(); orthologsIt++)
    {
      DNAIteratorConstPtr refDnaIt = *orthologsIt->first;
      set<DNAIteratorConstPtr *> *orthologSet = orthologsIt->second;
      char refDna = tolower(refDnaIt->getChar());
      if (refDna == 'n')
      {
        // Obviously shouldn't call snps here.
        continue;
      }
      hal_size_t numDifferentSpecies = 0; // # of species w/ base
                                          // different from ref
      for (set<DNAIteratorConstPtr *>::const_iterator orthologIt = orthologSet->begin(); orthologIt != orthologSet->end(); orthologIt++)
      {
        DNAIteratorConstPtr targetDnaIt = **orthologIt;
        char targetDna = tolower(targetDnaIt->getChar());
        if (targetDna == 'n')
        {
          continue;
        } else if (targetDna != refDna)
        {
          // This is a SNP for this species, but we have to wait until
          // the numDifferentSpecies is >= minSpeciesForSnp to call an
          // overall SNP.
          numDifferentSpecies++;
          (*numSnps)[targetDnaIt->getGenome()]++;
        }
        (*numOrthologousPairs)[targetDnaIt->getGenome()]++;
      }
//      printf("-----------\n");

      if (refTsvStream && numDifferentSpecies >= minSpeciesForSnp)
      {
        // Report a SNP to the TSV for this ortholog set.
        // First the sequence and position:
        const Sequence *refSeq = refDnaIt->getGenome()->getSequenceBySite(refDnaIt->getArrayIndex());
        refTsvStream << refSeq->getName() << "\t"
                     << refDnaIt->getArrayIndex() - refSeq->getStartPosition();
        // then the reference base:
        refTsvStream << "\t" << refDnaIt->getChar();
        // then finally the orthologs, in the same order that they
        // were spit out in the header.
        vector<char> orthologFields(targetGenomes.size());// initialized to '\0'
        for (set<DNAIteratorConstPtr *>::const_iterator i = orthologSet->begin();
             i != orthologSet->end(); i++)
        {
          char targetDna = (**i)->getChar();
          // - 3 because of the 3 fields preceding the target base fields
          assert(orthologFields[fieldForGenome[(**i)->getGenome()] - 3] == '\0');
          orthologFields[fieldForGenome[(**i)->getGenome()] - 3] = targetDna;
        }
        for (hal_size_t i = 0; i < orthologFields.size(); i++)
        {
          refTsvStream << "\t";
          if (orthologFields[i] != '\0')
          {
            refTsvStream << orthologFields[i];
          }
        }
        refTsvStream << endl;
      }
      delete orthologSet;
    }

    if (colIt->lastColumn()) {
      // Have to break here instead of at the beginning of the loop to
      // avoid missing the last column.
      break;
    }
    colIt->toRight();
  }
}
