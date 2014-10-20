// Replace ancestral nucleotides with the most likely given the tree
// and some substitution model. (assumes sites are independent)
#include "hal.h"
#include "sonLibTree.h"
#include "string.h"
#include "halBedScanner.h"
#include "ancestorsMLBed.h"
extern "C" {
#include "markov_matrix.h"
#include "tree_model.h"
}
using namespace std;
using namespace hal;

// globals -- for convenience
static bool writePosts; // for wigs
static double outValue; // for wigs

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal tree");
  optionsParser->addArgument("genome", "(ancestor) genome to modify");
  optionsParser->addArgument("model", "phyloP model file");
  optionsParser->addOption("writeHal", "write changes to hal file", false);
  optionsParser->addOption("startPos", "start position", 0);
  optionsParser->addOption("endPos", "end position", -1);
  optionsParser->addOption("sequence", "Sequence name. IMPORTANT: if "
                           "sequence name is not provided but startPos or "
                           "endPos are, they will be assumed to be in "
                           "genome coordinates", "");
  optionsParser->addOption("thresholdN", "threshold below which an N is output", 0.9);
  optionsParser->addOption("bed", "bed file to scan", "");
  optionsParser->addOptionFlag("outputPosts", "output posterior"
                               " probabilities for reference in wig"
                               " format", false);
  optionsParser->addOptionFlag("printWrites", "print base changes", false);
  return optionsParser;
}

typedef struct {
  const Genome *rootGenome;
  hal_index_t pos;
  // Reversed with respect to reference?
  bool reversed;
} rootInfo;

typedef struct {
  // Whether this node has already been calculated.
  bool done;
  // Position of this site in the genome
  hal_index_t pos;
  // Probability of leaves under this node given each nucleotide.
  double pLeaves[4];
  // This is a terrible and incorrect name
  // TODO change it without breaking everything
  double pOtherLeaves[4];
  // should only be set on the leaves at first.
  char dna;
  // Reversed with respect to reference?
  bool reversed;
  // phast ID from the model.
  int phastId;
  // Posterior probability of this call (in case we need it later)
  double post;
} felsensteinData;

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
          // TMP, unnecesssary but sanity check
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

/*
  double genRev(char childDNA, char parentDNA, MarkovMatrix *rateMatrix, double t)
  {
  MarkovMatrix *substMatrix;
  substMatrix = mm_new(4, "AGCT", DISCRETE);
  mm_exp(substMatrix, rateMatrix, t);
  return mm_get_by_state(substMatrix, parentDNA, childDNA);
  }*/
void doFelsenstein(stTree *node, TreeModel *mod)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(node);
  if (data->done) {
    return;
  }
  if (stTree_getChildNumber(node) == 0) {
    if (data->dna == 'N' || data->dna == 'n') {
      for (int dna = 0; dna < 4; ++dna) {
        data->pLeaves[dna] = 0.25;
      }
    } else {
      for (int dna = 0; dna < 4; ++dna) {
        if (dna == charToIndex(data->dna)) {
//          cout << stTree_getLabel(node) << ": " << data->dna << "( in " << (data->reversed ? "-" : "+") << " strand)" << endl;
          data->pLeaves[dna] = 1.0;
        } else {
          data->pLeaves[dna] = 0.0;
        }
      }
    }
  } else {
    for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
      doFelsenstein(stTree_getChild(node, childIdx), mod);
    }

    for (int dna = 0; dna < 4; dna++) {
      double prob = 1.0;
/*      // Need all possiblities for child nodes, so iterate over the
// cartesian product
// FIXME do this sensibly
assert(stTree_getChildNumber(node) < 16);
      
// get 4**stTree_getChildNumber(node)
int64_t cartesianLimit = 0;
for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
cartesianLimit *= 4;
}
for (int64_t cartesian = 0; cartesian < pow(4, stTree_getChildNumber(node)); cartesian++) {
int64_t c = cartesian;
double probAssignment = 1.0;
for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
stTree *childNode = stTree_getChild(node, childIdx);
felsensteinData *childData = (felsensteinData *) stTree_getClientData(childNode);
probAssignment *= probTransition(mod, childData->phastId, data->phastId, indexToChar(c % 4), indexToChar(dna));
probAssignment *= childData->pLeaves[c % 4];
c /= 4;
}
//        cout << "probability of assignment " << cartesian << ": " << probAssignment << endl;
prob += probAssignment;
}
*/
      for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
        // sum over the possibile assignments for this node
        stTree *childNode = stTree_getChild(node, childIdx);
        felsensteinData *childData = (felsensteinData *) stTree_getClientData(childNode);
        double probSubtree = 0.0;
        for (int childDna = 0; childDna < 4; childDna++) {
          probSubtree += childData->pLeaves[childDna]*probTransition(mod, childData->phastId, data->phastId, indexToChar(childDna), indexToChar(dna));
        }
        prob *= probSubtree;
      }
      data->pLeaves[dna] = prob;
//      cout << stTree_getLabel(node) << ": " << dna << " " << prob << endl;
      // Leave pointers to the letters in each child that are most
      // probable given each letter in this node

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
      double totalProb = 0.0;
      // trick found from phast code -- saves us some compute time
      double temp[4];
      for (int thisDna = 0; thisDna < 4; thisDna++) {
        temp[thisDna] = 0.0;
        for (int64_t j = 0; j < stTree_getChildNumber(tree); j++) {
          if (i == j) {
            continue;
          }
          stTree *siblingNode = stTree_getChild(tree, j);
          felsensteinData *siblingData = (felsensteinData *) stTree_getClientData(siblingNode);
          for (int siblingDna = 0; siblingDna < 4; siblingDna++) {
            temp[thisDna] += data->pOtherLeaves[thisDna] * siblingData->pLeaves[siblingDna] * probTransition(mod, siblingData->phastId, data->phastId, indexToChar(siblingDna), indexToChar(thisDna));
            //            cout << temp[thisDna] << " " << data->pOtherLeaves[thisDna] << " " << siblingData->pLeaves[siblingDna] << endl;
          }
        }
        if (stTree_getChildNumber(tree) == 1) {
          // Special case -- the sibling isn't in this tree.
          temp[thisDna] = data->pOtherLeaves[thisDna];
        }
      }
      for (int childDna = 0; childDna < 4; childDna++) {
        childData->pOtherLeaves[childDna] = 0.0;
        for (int thisDna = 0; thisDna < 4; thisDna++) {
          childData->pOtherLeaves[childDna] += temp[thisDna] * probTransition(mod, childData->phastId, data->phastId, indexToChar(childDna), indexToChar(thisDna));
        }
        totalProb += childData->pOtherLeaves[childDna] * childData->pLeaves[childDna];
      }
      int maxDna = -1;
      double maxProb = 0.0;
      for (int childDna = 0; childDna < 4; childDna++) {
        double post = childData->pOtherLeaves[childDna] * childData->pLeaves[childDna] / totalProb;
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
//      cout << string(1, childAssignment) << " : " << maxProb << endl;
      childData->post = maxProb;
      if (maxProb < threshold) {
//        cout << "probability " << maxProb << " below threshold for genome " << stTree_getLabel(childNode) << " at pos " << childData->pos << " original: " << string(1, childAssignment) << " different" << endl;
        childAssignment = 'N';
      }
      walkFelsenstein(mod, childNode, childAssignment, threshold);
    }
  }
}

void writeNucleotides(stTree *tree, AlignmentPtr alignment,
                      Genome *target, hal_index_t targetPos, bool writeHal = false,
                      bool printWrites = false)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  if (stTree_getChildNumber(tree) == 0) {
    return;
  }
  Genome *genome = alignment->openGenome(stTree_getLabel(tree));
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
    if (writeHal) {
      dnaIt->setChar(data->dna);
    }
  }
  if (writePosts && genome == target && data->pos == targetPos) {
    // correct genome and correct position
    outValue = data->post;
  }
  for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
    stTree *childNode = stTree_getChild(tree, i);
    writeNucleotides(childNode, alignment, target, targetPos, writeHal, printWrites);
  }
}

// TODO: just free the tree here as well, it'd be cleaner
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

void reEstimate(TreeModel *mod, AlignmentPtr alignment, Genome *genome, hal_index_t startPos, hal_index_t endPos, map<string, int> &nameToId, double threshold, bool writeHal, bool printWrites)
{
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
        outValue = 0.0;
        cout << outValue << endl;
      }
      continue;
    }
    doFelsenstein(tree, mod);
    // Find assignment for root node that maximizes P(leaves)
    felsensteinData *rootData = (felsensteinData *) stTree_getClientData(tree);
    // For prob(tree|char) -> prob(char|tree) (there is only one possible tree)
    double totalProbTree = 0.0;
    double maxProb = 0.0;
    int maxDna = -1;
    for (int dna = 0; dna < 4; dna++) {
      rootData->pOtherLeaves[dna] = 0.25;
      totalProbTree += rootData->pLeaves[dna];
      if (rootData->pLeaves[dna] > maxProb) {
        maxDna = dna;
        maxProb = rootData->pLeaves[dna];
      }
    }
    rootData->post = maxProb/totalProbTree;
    char assignment;
    if (maxDna == -1) {
      assignment = randNuc();
    } else if (rootData->post < threshold) {
//      cout << "below threshold for genome " << stTree_getLabel(tree) << " at pos " << rootData->pos << " original: " << string(1, indexToChar(maxDna)) << " different" << endl;
      assignment = 'N';
    } else {
      assignment = indexToChar(maxDna);
    }
    walkFelsenstein(mod, tree, assignment, threshold);
    writeNucleotides(tree, alignment, genome, pos, writeHal, printWrites);
    freeClientData(tree);
    stTree_destruct(tree);
    if (writePosts) {
      cout << outValue << endl;
    }
  }
}

int main(int argc, char *argv[])
{
  string halPath, genomeName, modPath, sequenceName, bedPath;
  CLParserPtr optParser = initParser();
  bool writeHal = false, printWrites = false;
  hal_index_t startPos = 0;
  hal_index_t endPos = -1;
  double threshold = 0.0;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    genomeName = optParser->getArgument<string>("genome");
    modPath = optParser->getArgument<string>("model");
    writeHal = optParser->getOption<bool>("writeHal");
    startPos = optParser->getOption<hal_index_t>("startPos");
    endPos = optParser->getOption<hal_index_t>("endPos");
    threshold = optParser->getOption<double>("thresholdN");
    sequenceName = optParser->getOption<string>("sequence");
    bedPath = optParser->getOption<string>("bed");
    writePosts = optParser->getFlag("outputPosts");
    printWrites = optParser->getFlag("printWrites");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }

  // Load phast model.
  FILE *infile = phast_fopen(modPath.c_str(), "r");
  TreeModel *mod = tm_new_from_file(infile, TRUE);
  phast_fclose(infile);
  tm_set_subst_matrices(mod);
  // Map names to phast model IDs.
  map<string, int> nameToId;
  List *phastList = tr_postorder(mod->tree);
  for (int i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = (TreeNode*) lst_get_ptr(phastList, i);
    nameToId[n->name] = n->id;
  }
  lst_free(phastList);

  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
    throw hal_exception("Genome " + genomeName + " not found in alignment.");
  }
  if (genome->getNumChildren() == 0) {
    throw hal_exception("Genome " + genomeName + " is a leaf genome.");
  }
  
  if (bedPath != "") {
    AncestorsMLBed bedScanner(mod, alignment, genome, nameToId, threshold,
                              writeHal, printWrites);
    bedScanner.scan(bedPath, -1);
    return 0;
  }

  if (sequenceName != "") {
    Sequence *sequence = genome->getSequence(sequenceName);
    if (sequence == NULL) {
      throw hal_exception("Sequence name not found!");
    }
    startPos += sequence->getStartPosition();
    if (endPos == -1) {
      endPos = sequence->getEndPosition();
    } else {
      endPos += sequence->getStartPosition();
      if (endPos > sequence->getEndPosition()) {
        endPos = sequence->getEndPosition();
      }
    }
  }

  if (endPos == -1 || endPos > genome->getSequenceLength()) {
    endPos = genome->getSequenceLength();
  }
  reEstimate(mod, alignment, genome, startPos, endPos, nameToId, threshold, writeHal, printWrites);
  alignment->close();
//  tm_free(mod);
  return 0;
}
