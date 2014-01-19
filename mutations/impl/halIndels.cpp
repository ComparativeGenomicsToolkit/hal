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
  // FIXME: pretty dumb way of asking for the normalization
  // denominator
  optionsParser->addOptionFlag("potentialSites",
                               "Output the total number of sites"
                               " that could be considered a clean indel if"
                               " there were an insertion/deletion there",
                               false);
  return optionsParser;
}

// check (inclusive) interval startPos--endPos for Ns.
static bool regionIsNotAmbiguous(const Genome* genome, hal_index_t startPos,
                                 hal_index_t endPos)
{
  if (startPos > endPos) {
    swap(startPos, endPos);
  }
  DNAIteratorConstPtr dnaIt = genome->getDNAIterator(startPos);
  hal_index_t length = endPos - startPos;
  for (hal_size_t i = 0; i < (hal_size_t) length; i++) {
    if (dnaIt->getChar() == 'N') {
      return false;
    }
    dnaIt->toRight();
  }
  return true;
}

static bool deletionIsNotAmbiguous(const ColumnIterator::ColumnMap *colMap,
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
    if (genome == refGenome) {
      continue;
    }
    const ColumnIterator::DNASet *dnaSet = colMapIt->second;
    assert(dnaSet->size() == 1);
    DNAIteratorConstPtr dnaIt = dnaSet->at(0);
    hal_index_t currPos = dnaIt->getArrayIndex();
    hal_index_t prevPos = (*prevPositions)[genome];
    (*prevPositions)[genome] = currPos;
    if (!regionIsNotAmbiguous(genome, currPos, prevPos)) {
      return false;
    }
  }
  return true;
}

// Check for Ns in any of the (strict single copy) targets
// FIXME: why are these named so that there end up being double- and
// triple-negatives in if conditionals
static bool isNotAmbiguous(const ColumnIterator::ColumnMap *colMap)
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

static void updatePrevPos(const ColumnIterator::ColumnMap *colMap,
                          map<const Genome *, hal_index_t> *prevPositions)
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
    (*prevPositions)[genome] = currPos;
  }
}

// returns true if the column is consistent with the previous positions
// Might crash if the column isn't strictly single copy.
static bool isContiguous(const ColumnIterator::ColumnMap *colMap,
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
static bool isStrictSingleCopy(const ColumnIterator::ColumnMap *colMap,
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

// FIXME: a lot of duplication from printIndels--but tough to separate
// out since there's a ton of insertion/deletion code we don't care
// about here
static void printNumCandidateSites(AlignmentConstPtr alignment,
                                   const Genome *refGenome,
                                   const set<const Genome *> targets,
                                   hal_size_t adjacentBases)
{
  hal_size_t refLength = refGenome->getSequenceLength();
  hal_size_t numSites = 0;
  ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets);
  PositionCache knownGoodSites;
  for (hal_index_t refPos = adjacentBases; refPos < refLength - adjacentBases;
       refPos++) {
    hal_index_t start = refPos - adjacentBases;
    hal_index_t end = refPos + adjacentBases;
    colIt->toSite(start, end, true);
    map <const Genome *, hal_index_t> prevPos;
    bool failedFiltering = false;
    while(1) {
      hal_index_t refColPos = colIt->getReferenceSequencePosition() + 
        colIt->getReferenceSequence()->getStartPosition();
      const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
      if (!knownGoodSites.find(refColPos)) {
        if (!isStrictSingleCopy(colMap, &targets) ||
            !isContiguous(colMap, &prevPos, 1, refGenome) ||
            !isNotAmbiguous(colMap)) {
          // we know this column doesn't pass filtering, so skip all
          // positions that we know will fail
          refPos = refColPos + adjacentBases; // 1 more will be added by for loop
          failedFiltering = true;
          break;
        } else {
          knownGoodSites.insert(refColPos);
        }
      }
      updatePrevPos(colMap, &prevPos);
      if (colIt->lastColumn()) {
        break;
      }
      colIt->toRight();
    }
    if (!failedFiltering) {
      numSites++;
    }
  }
  cout << numSites << endl;
}

static void printIndels(AlignmentConstPtr alignment, const Genome *refGenome,
                        const set<const Genome *> targets,
                        hal_size_t adjacentBases)
{
  RearrangementPtr rearrangement = refGenome->getRearrangement(0, 0, 1.0);
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
        const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
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
              if (!regionIsNotAmbiguous(refGenome, refGenomePos, breakEnd)) {
                // N in insertion. This could be a gap in a scaffold
                // so it's not considered clean.
                failedFiltering = true;
                break;
              }
              // Need to skip over the insertion.
              colIt->toSite(breakEnd + 1, end);
            }
          } else {
            break;
          }
        }
        if (!isStrictSingleCopy(colMap, &targets) ||
            !isContiguous(colMap, &prevPos, prevStep, refGenome) ||
            !isNotAmbiguous(colMap) ||
            // check that deleted regions don't have Ns (to keep
            // reversibility of insertion/deletion definitions)
            (prevStep != 1 && !deletionIsNotAmbiguous(colMap, &prevPos,
                                                      prevStep, refGenome))) {
          failedFiltering = true;
          break;
        }
        updatePrevPos(colMap, &prevPos);
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

int main(int argc, char *argv[])
{
  CLParserPtr optionsParser = initParser();
  string halPath, refGenomeName;
  hal_size_t adjacentBases;
  bool potentialSites;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    adjacentBases = optionsParser->getOption<hal_size_t>("adjacentBases");
    potentialSites = optionsParser->getFlag("potentialSites");
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

  // create target set: (ref, sibling(s), outgroup(s)).
  //
  // FIXME: In non-binary trees this will enforce the constraints on
  // *all* siblings, outgroups. This isn't what we want at
  // all--instead it should enforce the constraints on *at least* one
  // from each of siblings, outgroups. Otherwise deletions/insertions
  // shared in 2 of 3 children are not "clean".
  set <const Genome *> targets;
  const Genome *parentGenome = refGenome->getParent();
  for (hal_size_t i = 0; i < parentGenome->getNumChildren(); i++) {
    targets.insert(parentGenome->getChild(i));
  }
  if (refGenome->getParent()->getParent() != NULL) {
    // Add outgroup if possible (not child of root).
    const Genome *gpGenome = parentGenome->getParent();
    for (hal_size_t i = 0; i < parentGenome->getNumChildren(); i++) {
      if (gpGenome->getChild(i) != parentGenome) {
        targets.insert(gpGenome->getChild(i));
      }
    }
  }
  if (potentialSites) {
    printNumCandidateSites(alignment, refGenome, targets, adjacentBases);
  } else {
    printIndels(alignment, refGenome, targets, adjacentBases);
  }
 }
