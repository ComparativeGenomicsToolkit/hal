// Find clean indels
// TODO: merge into halBranchMutations
#include "hal.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance();
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("refGenome", 
                             "name of reference genome.");
  optionsParser->addOption("adjacentBases", "# of adjacent bases to examine "
                           "while filtering", 5);
  optionsParser->setDescription("Count (filtered) insertions/deletions in the "
                                "branch above the reference genome.");
  return optionsParser;
}

// Check for Ns in any of the (strict single copy) targets
bool isNotAmbiguous(const ColumnIterator::ColumnMap *colMap)
{
  ColumnIterator::ColumnMap::const_iterator colMapIt;
  for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++) {
    if (colMapIt->second->empty()) {
      // The column map can contain empty entries.
      continue;
    }
    const ColumnIterator::DNASet *dnaSet = colMapIt->second;
    assert(dnaSet->size() == 1);
    DNAIteratorConstPtr dnaIt = dnaSet->at(0);
    if (dnaIt->getChar() == 'N') {
      return false;
    }
  }
  return true;
}

// returns true if the column is consistent with the previous positions
// Might crash if the column isn't strictly single copy.
// Not obvious from the name, but a side effect is updating prevPos.
bool isContiguous(const ColumnIterator::ColumnMap *colMap,
                  map<const Genome *, hal_index_t> *prevPositions,
                  hal_size_t step, const Genome *refGenome)
{
  ColumnIterator::ColumnMap::const_iterator colMapIt;
  for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++) {
    if (colMapIt->second->empty()) {
      // The column map can contain empty entries.
      continue;
    }
    const Genome *genome = colMapIt->first->getGenome();
    const ColumnIterator::DNASet *dnaSet = colMapIt->second;
    assert(dnaSet->size() == 1);
    DNAIteratorConstPtr dnaIt = dnaSet->at(0);
    hal_index_t currPos = dnaIt->getArrayIndex();
    if (!prevPositions->count(genome)) {
      // initialize previous position map
      (*prevPositions)[genome] = currPos;
      continue;
    }
    hal_index_t prevPos = (*prevPositions)[genome];
    (*prevPositions)[genome] = currPos;
    // hacky. but the ref always steps by 1 even in deletions (obviously)
    hal_size_t myStep = (genome == refGenome) ? 1 : step;
    if (
      // Not adjacent in genome coordinates
      (dnaIt->getReversed() && currPos != prevPos - (hal_index_t) myStep) ||
      (!dnaIt->getReversed() && currPos != prevPos + (hal_index_t) myStep) ||
      // Not on same chromosome
      (genome->getSequenceBySite(currPos) != genome->getSequenceBySite(prevPos))
      ) {
      return false;
    }
  }
  return true;
}

// returns true if the column has exactly one entry for each genome.
//
// TODO: Should eventually merge w/ the crap in findSingleCopyRegions...
bool isStrictSingleCopy(const ColumnIterator::ColumnMap *colMap,
                        const set<const Genome *> *targets)
{
  ColumnIterator::ColumnMap::const_iterator colMapIt;
  set <const Genome *> seenGenomes;
  for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++) {
    if (colMapIt->second->empty()) {
      // The column map can contain empty entries.
      continue;
    }
    const Genome *colGenome = colMapIt->first->getGenome();
    if (seenGenomes.count(colGenome) || colMapIt->second->size() > 1) {
      return false;
    }
    seenGenomes.insert(colGenome);
  }
  if (seenGenomes.size() == targets->size()) {
#ifndef NDEBUG
    // Should be enough just to check the size. But we'll be careful--
    for (set <const Genome *>::iterator setIt = targets->begin();
         setIt != targets->end(); setIt++) {
      assert(seenGenomes.count(*setIt));
    }
#endif
    return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  CLParserPtr optionsParser = initParser();
  string halPath, refGenomeName;
  hal_size_t adjacentBases;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    adjacentBases = optionsParser->getOption<hal_size_t>("adjacentBases");
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser->printUsage(cerr);
    exit(1);
  }

  AlignmentConstPtr alignment = openHalAlignment(halPath, optionsParser);
  const Genome *refGenome = alignment->openGenome(refGenomeName);
  if (refGenome == NULL) {
    throw hal_exception("Genome " + refGenomeName + " does not exist");
  }
  if (refGenome->getParent() == NULL) {
    throw hal_exception("Cannot use the root genome as a reference.");
  }
  RearrangementPtr rearrangement = refGenome->getRearrangement(0, 0, 1.0);

  // create target set: (ref, sibling(s), outgroup(s)).
  //
  // FIXME: In non-binary trees this will enforce the constraints on
  // *all* siblings, outgroups. This isn't what we want at
  // all--instead it should enforce the constraints on *at least* one
  // from each of siblings, outgroups. Otherwise deletions/insertions
  // shared in 2 of 3 children are not "clean".
  set <const Genome *> targets;
  if (refGenome->getParent() != NULL &&
      refGenome->getParent()->getParent() != NULL) {
    const Genome *parentGenome = refGenome->getParent();
    for (hal_size_t i = 0; i < parentGenome->getNumChildren(); i++) {
      targets.insert(parentGenome->getChild(i));
    }
    const Genome *gpGenome = parentGenome->getParent();
    for (hal_size_t i = 0; i < parentGenome->getNumChildren(); i++) {
      if (gpGenome->getChild(i) != parentGenome) {
        targets.insert(gpGenome->getChild(i));
      }
    }
  }

  do {
    if (rearrangement->getID() == Rearrangement::Insertion ||
        rearrangement->getID() == Rearrangement::Deletion) {
      hal_index_t breakStart = rearrangement->getLeftBreakpoint()->getStartPosition();
      hal_index_t breakEnd = rearrangement->getRightBreakpoint()->getEndPosition();
      if (rearrangement->getID() == Rearrangement::Deletion) {
        // The right breakpoint seems to not be set or something in deletions.
        // Since the left breakpoint = right breakpoint in deletions,
        // set it manually
        // TODO: fix rearrangement class
        breakEnd = breakStart;
      }
      hal_index_t start = breakStart - adjacentBases;
      hal_index_t end = breakEnd + adjacentBases;
      if (start < 0 || end >= (hal_index_t) refGenome->getSequenceLength()) {
        // Indels very close to the end of sequences can't be "clean."
        continue;
      }
      ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets, 0,
                                                             start, end);
      map <const Genome *, hal_index_t> prevPos;
      bool failedFiltering = false;
      hal_size_t prevStep = 32432432423; // just to catch bugs...
      while(1) {
        hal_index_t refGenomePos = colIt->getReferenceSequencePosition() + 
          colIt->getReferenceSequence()->getStartPosition();
        if (refGenomePos == breakStart) {
          // Fiddle with the prevPos map so we only enforce the
          // adjacencies we need (adjacency in the reference for a
          // deletion, adjacencies in all other genomes for an
          // insertion.)
          if (breakEnd + 1 < end) {
            if (rearrangement->getID() == Rearrangement::Deletion) {
              rearrangement->identifyDeletionFromLeftBreakpoint(rearrangement->getLeftBreakpoint());
              pair<hal_index_t, hal_index_t> deletedRange = rearrangement->getDeletedRange();
              prevStep = llabs(deletedRange.first - deletedRange.second) + 2;
            } else {
              // just in case the condition is changed above
              assert(rearrangement->getID() == Rearrangement::Insertion);
              prevPos.erase(refGenome);
              // Need to skip over the insertion.
              colIt->toSite(breakEnd + 1, end);
            }
          } else {
            break;
          }
        }
        const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();

        if (!isStrictSingleCopy(colMap, &targets) ||
            !isContiguous(colMap, &prevPos, prevStep, refGenome) ||
            !isNotAmbiguous(colMap)) {
          failedFiltering = true;
          break;
        }
        if (colIt->lastColumn()) {
          break;
        }
        colIt->toRight();
        prevStep = 1;
      }
      if (!failedFiltering) {
        if (rearrangement->getID() == Rearrangement::Deletion) {
          pair<hal_index_t, hal_index_t> deletedRange = rearrangement->getDeletedRange();
          const Genome *parent = refGenome->getParent();
          const Sequence *seq = parent->getSequenceBySite(deletedRange.first);
          assert(seq == parent->getSequenceBySite(deletedRange.second));
          cout << seq->getName() << "\t"
               << deletedRange.first - seq->getStartPosition() << "\t"
               << deletedRange.second - seq->getStartPosition() + 1 << "\tD\t"
               << parent->getName() << "\t" << refGenome->getName() << endl;
        } else {
          const Genome *parent = refGenome->getParent();
          const Sequence *seq = refGenome->getSequenceBySite(start);
          assert(seq == refGenome->getSequenceBySite(end));
          cout << seq->getName() << "\t"
               << breakStart - seq->getStartPosition() << "\t"
               << breakEnd - seq->getStartPosition() + 1 << "\tI\t"
               << parent->getName() << "\t" << refGenome->getName() << endl;
        }
      }
    }
  } while(rearrangement->identifyNext() == true);
}
