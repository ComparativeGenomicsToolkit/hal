// Find clean indels
// TODO: merge into halBranchMutations
#include "hal.h"

using namespace std;
using namespace hal;

// TODO: sloppy, but whatever -- this is all going to be rearranged to
// fit in halBranchMutations soon anyway
enum indelType {
  NONE,
  INSERTION,
  DELETION
};

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
  optionsParser->addOptionFlag("onlyExtantTargets",
                               "Use only extant genomes for 'sibling'/outgroup",
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
                                   const Genome *refGenome)
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

// report if this is a (potentially unclean) insertion in the
// reference relative to the other targets
static bool isInsertion(const ColumnIterator::ColumnMap *colMap,
                        const Genome *refGenome)
{
  ColumnIterator::ColumnMap::const_iterator colMapIt;
  hal_size_t numCopies = 0;
  for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++) {
    if (colMapIt->second->empty()) {
      // The column map can contain empty entries.
      continue;
    }
    const Genome *colGenome = colMapIt->first->getGenome();
    if (colGenome != refGenome) {
      // since we are only traversing the targets this is OK to do, if
      // we are traversing ancestors or something like that it could
      // be problematic
      return false;
    }
    numCopies++;
  }
  if (numCopies > 1) {
    return false;
  }
  return true;
}

// report deletion size if this is a (potentially unclean) deletion
// relative to the other targets
// otherwise 0
static hal_size_t getDeletedSize(const ColumnIterator::ColumnMap *colMap,
                                 map<const Genome *, hal_index_t> *prevPositions,
                                 const Genome *refGenome)
{
  ColumnIterator::ColumnMap::const_iterator colMapIt;
  hal_size_t delSize = 0;
  for (colMapIt = colMap->begin(); colMapIt != colMap->end(); colMapIt++) {
    if (colMapIt->second->empty()) {
      // The column map can contain empty entries.
      continue;
    }
    const Genome *colGenome = colMapIt->first->getGenome();
    const ColumnIterator::DNASet *dnaSet = colMapIt->second;
    assert(dnaSet->size() == 1);
    DNAIteratorConstPtr dnaIt = dnaSet->at(0);
    hal_index_t currPos = dnaIt->getArrayIndex();
    if (!prevPositions->count(colGenome)) {
      // initialize previous position map
      (*prevPositions)[colGenome] = currPos;
      continue;
    }
    hal_index_t prevPos = (*prevPositions)[colGenome];
    if (colGenome == refGenome) {
      assert(currPos = prevPos + 1);
    }
    if (
      ((dnaIt->getReversed() && currPos != prevPos - 1) ||
       (!dnaIt->getReversed() && currPos != prevPos + 1)) &&
      (colGenome->getSequenceBySite(currPos) == colGenome->getSequenceBySite(prevPos))
      ) {
      if (delSize != 0) {
        // There has already been a deletion in another target, check
        // that they are the same length
        hal_size_t myDelSize = llabs(currPos - prevPos) - 1;
        assert(myDelSize > 0);
        if (delSize != myDelSize) {
          // Disagreement on deletion length between sister &
          // outgroup, so this can never be a clean deletion
          return 0;
        }
      } else {
        // initialize deletion length -- -1 because currPos is 1 past
        // the deletion
        delSize = llabs(currPos - prevPos) - 1;
        assert(delSize > 0);
      }
    } else {
      // Not deleted
      if (colGenome != refGenome) {
        // Not deleted in all the other targets, so for our purposes
        // not deleted at all.
        return 0;
      }
    }
  }
  return delSize;
}

// get information about an indel, which starts at refPos, if one is present.
static pair<indelType, hal_size_t>
getIndel(hal_index_t refPos,
         const Genome *refGenome,
         const set<const Genome *> *targets)
{
  if (refPos == 0) {
    return make_pair(NONE, 0);
  }
  ColumnIteratorPtr colIt = refGenome->getColumnIterator(targets, 0,
                                                         refPos - 1);
  const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
  if (!isStrictSingleCopy(colMap, targets)) {
    // Make sure our assumptions hold about prevPos maps
    return make_pair(NONE, 0);
  }
  map<const Genome *, hal_index_t> prevPos;
  updatePrevPos(colMap, &prevPos);
  colIt->toRight();
  colMap = colIt->getColumnMap();
  // if current base is not present in the other targets eat up sequence
  // until end of insertion, call unclean insertion of length X
  if (isInsertion(colMap, refGenome)) {
    while (isInsertion(colMap, refGenome)) {
      colIt->toRight();
      colMap = colIt->getColumnMap();
      if (colIt->lastColumn()) {
        // impossible to call clean insertion at end of genome.
        return make_pair(NONE, 0);
      }
    }
    hal_index_t currPos = colIt->getReferenceSequencePosition()
      + colIt->getReferenceSequence()->getStartPosition();
    assert(currPos > refPos);
    hal_size_t insertedSize = currPos - refPos;
    if (!regionIsNotAmbiguous(refGenome, refPos, refPos + insertedSize)) {
      // N in insertion. This could be a gap in a scaffold so it's not
      // considered clean.
      return make_pair(NONE, 0);
    }
    return make_pair(INSERTION, insertedSize);
  }

  // if this base skips X bases in both the other targets call an
  // unclean deletion of length X
  if (!isStrictSingleCopy(colMap, targets)) {
    // Make sure our assumptions in getDeletedSize hold for this
    // column
    return make_pair(NONE, 0);
  }
  hal_size_t deletedSize = getDeletedSize(colMap, &prevPos, refGenome);
  if (deletedSize) {
    return make_pair(DELETION, deletedSize);
  }

  return make_pair(NONE, 0);
}

static void printIndels(const Genome *refGenome,
                        const set<const Genome *> targets,
                        hal_size_t adjacentBases)
{
  hal_size_t refLength = refGenome->getSequenceLength();
  hal_size_t numSites = 0;
  ColumnIteratorPtr colIt = refGenome->getColumnIterator(&targets);
  // good flanking site
  PositionCache knownGoodSites;
  for (hal_index_t refPos = adjacentBases; refPos < refLength - adjacentBases;
       refPos++) {
    pair<indelType, hal_size_t> indel;
    indel = getIndel(refPos, refGenome, &targets);
    hal_index_t start = refPos - adjacentBases;
    hal_index_t end = refPos + adjacentBases;
    if (indel.first == INSERTION) {
      end += indel.second;
    }
    colIt->toSite(start, end, true);
    map <const Genome *, hal_index_t> prevPos;
    bool failedFiltering = false;
    hal_size_t step = 1;
    while(1) {
      hal_index_t refColPos = colIt->getReferenceSequencePosition() + 
        colIt->getReferenceSequence()->getStartPosition();
      if (refColPos == refPos && indel.first == DELETION) {
        // jump "step" bases -- i.e. past the deleted region
        step = indel.second + 1;
      } else if (refColPos == refPos && indel.first == INSERTION) {
        // don't enforce adjacency on insertion since we're skipping it
        prevPos.erase(refGenome);
        if (refGenome->getSequenceBySite(refPos) !=
            refGenome->getSequenceBySite(refPos + indel.second)) {
          // Insertion crosses sequence end
          failedFiltering = true;
          break;
        }
        colIt->toSite(refPos + indel.second, end);
        continue;
      } else {
        step = 1;
      }
      const ColumnIterator::ColumnMap *colMap = colIt->getColumnMap();
      if (!knownGoodSites.find(refColPos)) {
        if (!isStrictSingleCopy(colMap, &targets) ||
            !isContiguous(colMap, &prevPos, step, refGenome) ||
            !isNotAmbiguous(colMap) ||
            (step != 1 && !deletionIsNotAmbiguous(colMap, &prevPos,
                                                  refGenome))) {
          failedFiltering = true;
          if (indel.first == INSERTION) {
            // failed indel means that we don't have to check anywhere
            // inside the insertion, it will automatically fail
            refPos += indel.second;
          }
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
    if (indel.first != NONE && !failedFiltering) {
      if (indel.first == DELETION) {
        const Sequence *seq = refGenome->getSequenceBySite(refPos);
        cout << seq->getName() << "\t"
             << refPos - seq->getStartPosition() << "\t"
             << refPos - seq->getStartPosition()<< "\tD\t" 
             << indel.second << endl;
      } else {
        const Sequence *seq = refGenome->getSequenceBySite(refPos);
        assert(seq == refGenome->getSequenceBySite(refPos + indel.second));
        cout << seq->getName() << "\t"
             << refPos - seq->getStartPosition() << "\t"
             << refPos + indel.second - seq->getStartPosition() << "\tI\t"
             << endl;
        refPos += indel.second;
      }
    }
    if (!failedFiltering) {
      numSites++;
    }
  }
  cout << "# num sites possible: " << numSites << endl;
}

static pair<double, const Genome *>
findMinPathUnder(const Genome *parent,
                 const Genome *exclude,
                 double branchLengthSoFar = 0)
{
  if (parent->getNumChildren() == 0) {
    return make_pair(branchLengthSoFar, parent);
  }
  double bestBranchLength = INFINITY;
  const Genome *bestGenome = NULL;
  for (hal_size_t i = 0; i < parent->getNumChildren(); i++) {
    const Genome *child = parent->getChild(i);
    if (child == exclude) {
      continue;
    }
    const Alignment *alignment = parent->getAlignment();
    double branchLength = alignment->getBranchLength(parent->getName(),
                                                     child->getName());
    pair<double, const Genome *> bestLocalPath = findMinPathUnder(child,
                                                                  exclude,
                                                                  branchLengthSoFar + branchLength);
    if (bestLocalPath.first < bestBranchLength) {
      bestGenome = bestLocalPath.second;
      bestBranchLength = bestLocalPath.first;
    }
  }
  return make_pair(bestBranchLength, bestGenome);
}

int main(int argc, char *argv[])
{
  CLParserPtr optionsParser = initParser();
  string halPath, refGenomeName;
  hal_size_t adjacentBases;
  bool onlyExtantTargets;
  try
  {
    optionsParser->parseOptions(argc, argv);
    halPath = optionsParser->getArgument<string>("halFile");
    refGenomeName = optionsParser->getArgument<string>("refGenome");
    adjacentBases = optionsParser->getOption<hal_size_t>("adjacentBases");
    onlyExtantTargets = optionsParser->getFlag("onlyExtantTargets");
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

  set <const Genome *> targets;
  // create target set: (ref, sibling(s), outgroup(s)).
  if (onlyExtantTargets) {
    // Potentially our ref could have been used to influence insertion/
    // deletion calls in ancestral nodes. So we can use exclusively extant
    // genomes to avoid that.
    // FIXME: In non-binary trees it might be desirable to require agreement
    // among *any* outgroup or *any* "sibling" rather than just the closest
    // one.
    const Genome *parentGenome = refGenome->getParent();
    // Add closest "sibling" (actually closest extant genome under sibling)
    targets.insert(findMinPathUnder(parentGenome, refGenome).second);
    if (refGenome->getParent()->getParent() != NULL) {
      const Genome *gpGenome = parentGenome->getParent();
      // add closest outgroup.
      targets.insert(findMinPathUnder(gpGenome, parentGenome).second);
    }
    // Add reference (used elsewhere to check single-copy quickly)
    targets.insert(refGenome);
  } else {
    // FIXME: In non-binary trees this will enforce the constraints on
    // *all* siblings, outgroups. This isn't what we want at
    // all--instead it should enforce the constraints on *at least*
    // one, or maybe at least 2,3,...  from each of siblings,
    // outgroups. Otherwise deletions/insertions shared in 2 of 3
    // children are not "clean".
    const Genome *parentGenome = refGenome->getParent();
    // add siblings (and reference! This is needed elsewhere)
    for (hal_size_t i = 0; i < parentGenome->getNumChildren(); i++) {
      targets.insert(parentGenome->getChild(i));
    }
    if (refGenome->getParent()->getParent() != NULL) {
      // Add outgroup if possible (not child of root).
      const Genome *gpGenome = parentGenome->getParent();
      for (hal_size_t i = 0; i < gpGenome->getNumChildren(); i++) {
        if (gpGenome->getChild(i) != parentGenome) {
          targets.insert(gpGenome->getChild(i));
        }
      }
    }
  }
  printIndels(refGenome, targets, adjacentBases);
}
