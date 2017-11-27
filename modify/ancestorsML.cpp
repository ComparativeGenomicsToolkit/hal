// Replace ancestral nucleotides with the most likely given the tree
// and some substitution model. (assumes sites are independent)
#include "hal.h"
#include "sonLibTree.h"
#include "string.h"
#include "halBedScanner.h"
#include "ancestorsML.h"
extern "C" {
#include "markov_matrix.h"
#include "tree_model.h"
}
// PHAST code defines min, max macros which conflict with the reserved C++ names.
#undef min
#undef max
using namespace std;
using namespace hal;

// globals -- for convenience
static double outValue; // for wigs

// sum log-transformed probabilities.
static inline double log_space_add(double x, double y) {
  if (x == -INFINITY) {
    return y;
  }
  if (y == -INFINITY) {
    return x;
  }
  return max(x, y) + log1p(exp(-fabs(x - y)));
}

// These could both be replaced with global arrays
inline char indexToChar(int index) {
  switch(index) {
  case 0: return 'A';
  case 1: return 'G';
  case 2: return 'C';
  case 3: return 'T';
  default: throw hal_exception("Invalid index");
  }
}

char randNuc(void)
{
  static char nucs[] = {'A', 'C', 'G', 'T'};
  return nucs[random() % 4];
}

inline int charToIndex(char dna)
{
  switch(dna) {
  case 'A':
  case 'a':
    return 0;
  case 'G':
  case 'g':
    return 1;
  case 'C':
  case 'c':
    return 2;
  case 'T':
  case 't':
    return 3;
  case 'N':
  case 'n':
    return -1;
  default:
    throw hal_exception("Unsupported character " + string(1, dna) + 
                        " found in sequence");
  }
}

// find root genome of tree for this position in the genome
rootInfo * findRoot(const Genome *genome,
                    hal_index_t pos,
                    bool reversed=false) {
  if (genome->getParent() == NULL) {
    rootInfo *data = (rootInfo *) malloc(sizeof(rootInfo));
    data->rootGenome = genome;
    data->pos = pos;
    data->reversed = reversed;
    return data;
  }
  TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
  topIt->toSite(pos, false);
  bool parentReversed = topIt->getParentReversed() ? !reversed : reversed;
  if (!topIt->hasParent()) {
//    cout << "Root genome for pos: " << pos << " " << genome->getName() << endl;
    rootInfo *data = (rootInfo *) malloc(sizeof(rootInfo));
    data->rootGenome = genome;
    data->pos = pos;
    data->reversed = reversed;
    return data;
  } else {
    BottomSegmentIteratorConstPtr botIt = genome->getParent()->getBottomSegmentIterator();
    botIt->toParent(topIt);
    hal_index_t parentPos = botIt->getStartPosition();
    hal_index_t offset = abs(pos - topIt->getStartPosition());
    parentPos = topIt->getParentReversed() ? parentPos - offset : parentPos + offset;
    //  cout << "Going to parent " << genome->getParent()->getName() << ": start pos " << pos << " parent pos: " << parentPos << "reversed: " << parentReversed << endl;
    return findRoot(genome->getParent(), parentPos, parentReversed);
  }
}

// Remove a node from the tree and prune any unneeded nodes left over
// (ancestral nodes that would have no children).
void removeAndPrune(stTree *tree) {
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  free(data);
  stTree *parentNode = stTree_getParent(tree);
  if (parentNode != NULL) {
    stTree_setParent(tree, NULL);
    // We need to avoid freeing the root so main can know that there
    // is no point in continuing for this site (root has no children)
    if (stTree_getChildNumber(parentNode) == 0
        && stTree_getParent(parentNode) != NULL) {
//      removeAndPrune(parentNode);
    }
  }
  stTree_destruct(tree);
}

// Build site-specific tree below this genome.
void buildTree(AlignmentConstPtr alignment, const Genome *genome,
               hal_index_t pos, stTree *tree, bool reversed, map<string, int> *nameToId = NULL)
{
  stTree_setLabel(tree, genome->getName().c_str());
  felsensteinData *data = (felsensteinData *) malloc(sizeof(felsensteinData));
  memset(data, 0, sizeof(felsensteinData));
  data->pos = pos;
  data->reversed = reversed;
  if (nameToId != NULL) {
    data->phastId = (*nameToId)[genome->getName()];
  }
  stTree_setClientData(tree, data);
  if (genome->getNumChildren() == 0) {
    DNAIteratorConstPtr dnaIt = genome->getDNAIterator(pos);
    if (reversed) {
      dnaIt->toReverse();
    }
    data->dna = dnaIt->getChar();
    return;
  }
  data->dna = 'Z'; // signals an ancestor for the pruning process -- hacky
  BottomSegmentIteratorConstPtr botIt = genome->getBottomSegmentIterator();
  botIt->toSite(pos, false);
  assert(botIt->getReversed() == false);
  for (hal_size_t i = 0; i < botIt->getNumChildren(); i++) {
    hal_index_t childIndex = botIt->getChildIndex(i);
    if (childIndex != NULL_INDEX) {
      const Genome *childGenome = genome->getChild(i);
      TopSegmentIteratorConstPtr topIt = childGenome->getTopSegmentIterator();
      topIt->toChild(botIt, i);
      if (topIt->getNextParalogyIndex() != NULL_INDEX) {
        // Go through the paralogy cycle and add the
        // paralogous sites.
        // NOTE!: can theoretically run into problems
        // comparing these iterators if there is an
        // orientation change between paralogous segments (can
        // possibly create more duplications than really
        // exist) Won't happen with the way the iterator
        // comparison is implemented now though.
        TopSegmentIteratorConstPtr original = topIt->copy();
        for (topIt->toNextParalogy(); !topIt->equals(original); topIt->toNextParalogy()) {
          // sanity check
          assert(topIt->getLength() == botIt->getLength());
          hal_index_t childPos = topIt->getStartPosition();
          hal_index_t offset = abs(pos - botIt->getStartPosition());
          childPos = topIt->getParentReversed() ? childPos - offset : childPos + offset;
          stTree *childNode = stTree_construct();
          stTree_setParent(childNode, tree);
          double branchLength = alignment->getBranchLength(genome->getName(),
                                                           childGenome->getName());
          // TODO: make option for branch length of duplications
          // (e.g. as fraction of species-tree branch length)
          stTree_setBranchLength(childNode, branchLength);
          bool childReversed = topIt->getParentReversed() ? !reversed : reversed;
          buildTree(alignment, childGenome, childPos, childNode, childReversed, nameToId);
        }
      }
      hal_index_t childPos = topIt->getStartPosition();
      hal_index_t offset = abs(pos - botIt->getStartPosition());
      childPos = botIt->getChildReversed(i) ? childPos - offset : childPos + offset;
      stTree *childNode = stTree_construct();
      stTree_setParent(childNode, tree);
      double branchLength = alignment->getBranchLength(genome->getName(),
                                                       childGenome->getName());
      stTree_setBranchLength(childNode, branchLength);
      bool childReversed = botIt->getChildReversed(i) ? !reversed : reversed;
      buildTree(alignment, childGenome, childPos, childNode, childReversed, nameToId);
    }
  }
}

bool pruneTree(stTree *tree)
{
  // This function removes any ancestral leaf (usually from alignment
  // slop aligning to the edge of an ancestral scaffold gap).
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
    if(pruneTree(stTree_getChild(tree, i))) {
      // The child numbers will shift by 1 on deletion.
      i -= 1;
    }
  }

  if (stTree_getChildNumber(tree) == 0 && stTree_getParent(tree) != NULL) {
    felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
    // hack, but check if this is an ancestor in the species tree.
    if (data->dna == 'Z') {
      removeAndPrune(tree);
      return true;
    }
  }
  return false;
}

double probTransition(TreeModel *mod, int childId, int parentId, char childDNA, char parentDNA)
{
  assert(mod->nratecats == 1);
  assert(mod->P[childId][0] != NULL);
  MarkovMatrix *substMatrix = mod->P[childId][0];
  return mm_get_by_state(substMatrix, parentDNA, childDNA);
}

void doFelsenstein(stTree *node, TreeModel *mod)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(node);
  if (data->done) {
    return;
  }
  if (stTree_getChildNumber(node) == 0) {
    if (data->dna == 'N' || data->dna == 'n') {
      for (int dna = 0; dna < 4; ++dna) {
        data->pLeaves[dna] = log(0.25);
      }
    } else {
      for (int dna = 0; dna < 4; ++dna) {
        if (dna == charToIndex(data->dna)) {
          data->pLeaves[dna] = log(1.0);
        } else {
          data->pLeaves[dna] = -INFINITY;
        }
      }
    }
  } else {
    for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
      doFelsenstein(stTree_getChild(node, childIdx), mod);
    }

    for (int dna = 0; dna < 4; dna++) {
      double prob = 0.0;
      for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
        // sum over the possibile assignments for this node
        stTree *childNode = stTree_getChild(node, childIdx);
        felsensteinData *childData = (felsensteinData *) stTree_getClientData(childNode);
        double probSubtree = -INFINITY;
        for (int childDna = 0; childDna < 4; childDna++) {
          double probBranch = log(probTransition(mod, childData->phastId, data->phastId,
                                                 indexToChar(childDna), indexToChar(dna)));
          probSubtree = log_space_add(probSubtree, childData->pLeaves[childDna] + probBranch);
        }
        prob += probSubtree;
      }
      data->pLeaves[dna] = prob;
    }
  }
  data->done = true;
}

// Assign nucleotides to each node in the tree.
void walkFelsenstein(TreeModel *mod, stTree *tree, char assignment, double threshold)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  data->dna = assignment;
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
    stTree *childNode = stTree_getChild(tree, i);
    if (stTree_getChildNumber(childNode) != 0) {
      felsensteinData *childData = (felsensteinData *) stTree_getClientData(childNode);

      // Find posterior probability of this base.
      double totalProb = -INFINITY;
      // trick found from phast code -- saves us some compute time
      double temp[4];
      for (int thisDna = 0; thisDna < 4; thisDna++) {
        temp[thisDna] = -INFINITY;
        for (int64_t j = 0; j < stTree_getChildNumber(tree); j++) {
          if (i == j) {
            continue;
          }
          stTree *siblingNode = stTree_getChild(tree, j);
          felsensteinData *siblingData = (felsensteinData *) stTree_getClientData(siblingNode);
          for (int siblingDna = 0; siblingDna < 4; siblingDna++) {
            temp[thisDna] = log_space_add(temp[thisDna], data->pOtherLeaves[thisDna] + siblingData->pLeaves[siblingDna] + log(probTransition(mod, siblingData->phastId, data->phastId, indexToChar(siblingDna), indexToChar(thisDna))));
          }
        }
        if (stTree_getChildNumber(tree) == 1) {
          // Special case -- the sibling isn't in this tree.
          temp[thisDna] = data->pOtherLeaves[thisDna];
        }
      }
      for (int childDna = 0; childDna < 4; childDna++) {
        childData->pOtherLeaves[childDna] = -INFINITY;
        for (int thisDna = 0; thisDna < 4; thisDna++) {
          childData->pOtherLeaves[childDna] = log_space_add(childData->pOtherLeaves[childDna], temp[thisDna] + log(probTransition(mod, childData->phastId, data->phastId, indexToChar(childDna), indexToChar(thisDna))));
        }
        totalProb = log_space_add(totalProb, childData->pOtherLeaves[childDna] + childData->pLeaves[childDna]);
      }
      int maxDna = -1;
      double maxProb = -INFINITY;
      for (int childDna = 0; childDna < 4; childDna++) {
        double post = childData->pOtherLeaves[childDna] + childData->pLeaves[childDna] - totalProb;
        if (post > maxProb) {
          maxDna = childDna;
          maxProb = post;
        }
      }
      char childAssignment;
      if (maxDna == -1) {
        childAssignment = randNuc();
      } else {
        childAssignment = indexToChar(maxDna);
      }
      childData->post = maxProb;
      if (maxProb < threshold) {
        childAssignment = 'N';
      }
      walkFelsenstein(mod, childNode, childAssignment, threshold);
    }
  }
}

void writeNucleotides(stTree *tree, AlignmentConstPtr alignment,
                      const Genome *target, hal_index_t targetPos,
                      bool printWrites = false)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  if (stTree_getChildNumber(tree) == 0) {
    return;
  }
  const Genome *genome = alignment->openGenome(stTree_getLabel(tree));
  assert(genome != NULL);
  DNAIteratorPtr dnaIt = genome->getDNAIterator(data->pos);
  if (data->reversed) {
    dnaIt->toReverse();
  }
  char dna = toupper(dnaIt->getChar());
  if (data->dna != dna) {
    if (printWrites) {
      cout << genome->getName() << "\t" << data->pos << "\t" << string(1, dna) << "\t" << string(1, data->dna) << endl;
    }
  }
  if (genome == target && data->pos == targetPos) {
    // correct genome and correct position
    outValue = data->post;
  }
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
    stTree *childNode = stTree_getChild(tree, i);
    writeNucleotides(childNode, alignment, target, targetPos, printWrites);
  }
}

void freeClientData(stTree *tree)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  if (stTree_getChildNumber(tree)) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
//        free(data->childrenCharProbs[i][j]);
      }
    }
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
      freeClientData(stTree_getChild(tree, i));
    }
  }
  free(data);
}

void reEstimate(TreeModel *mod, AlignmentConstPtr alignment, const Genome *genome, hal_index_t startPos, hal_index_t endPos, map<string, int> &nameToId, double threshold, bool printWrites, bool writePosts)
{
  threshold = log(threshold);
  stTree *tree = NULL;
  bool firstRun = true;
  for (hal_index_t pos = startPos; pos < endPos; pos++) {
    outValue = 0.0;
    if (writePosts && firstRun) {
      const Sequence *seq = genome->getSequenceBySite(pos);
      // position + 1 because wigs are 1-based.
      cout << "fixedStep chrom=" << seq->getName() << 
        " start=" << pos - seq->getStartPosition() + 1 << " step=1" << endl;
      firstRun = false;
    }
    tree = stTree_construct();
    // Find root of tree
    rootInfo *rootInfo = findRoot(genome, pos);
    const Genome *root = rootInfo->rootGenome;
    hal_index_t rootPos = rootInfo->pos;
    bool rootReversed = rootInfo->reversed;
    buildTree(alignment, root, rootPos, tree, rootReversed, &nameToId);
    pruneTree(tree);
    if (stTree_getChildNumber(tree) == 0) {
      // No reason to build a tree, there's an insertion in the root
      // node relative to its children.
      freeClientData(tree);
      stTree_destruct(tree);
      if (writePosts) {
        // need to keep the wig in order
        outValue = -INFINITY;
        cout << outValue << endl;
      }
      continue;
    }
    doFelsenstein(tree, mod);
    // Find assignment for root node that maximizes P(leaves)
    felsensteinData *rootData = (felsensteinData *) stTree_getClientData(tree);
    // For prob(tree|char) -> prob(char|tree) (there is only one possible tree)
    double totalProbTree = -INFINITY;
    double maxProb = -INFINITY;
    int maxDna = -1;
    for (int dna = 0; dna < 4; dna++) {
      rootData->pOtherLeaves[dna] = log(0.25);
      totalProbTree = log_space_add(totalProbTree, rootData->pLeaves[dna]);
      if (rootData->pLeaves[dna] > maxProb) {
        maxDna = dna;
        maxProb = rootData->pLeaves[dna];
      }
    }
    rootData->post = maxProb - totalProbTree;
    char assignment;
    if (maxDna == -1) {
      assignment = randNuc();
    } else if (rootData->post < threshold) {
      assignment = 'N';
    } else {
      assignment = indexToChar(maxDna);
    }
    walkFelsenstein(mod, tree, assignment, threshold);
    writeNucleotides(tree, alignment, genome, pos, printWrites);
    freeClientData(tree);
    stTree_destruct(tree);
    if (writePosts) {
      cout << outValue << endl;
    }
  }
}
