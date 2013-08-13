// Replace ancestral nucleotides with the most likely given the tree
// and some substitution model. (assumes sites are independent)
//
// Currently broken:
//  - Doesn't handle ancestors with insertions with respect to all their
//    children.
#include "hal.h"
#include "sonLibTree.h"
#include "string.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal tree");
  optionsParser->addArgument("genome", "(ancestor) genome to modify");
  optionsParser->addArgument("substRate", "jukes-cantor substitution rate");
  return optionsParser;
}

typedef struct {
  // Whether this node has already been calculated.
  bool done;
  // Position of this site in the genome
  hal_index_t pos;
  // Probability of leaves under this node given each nucleotide.
  double pLeaves[4];
  // points to maximum probability chars for the children of this node,
  // given each nucleotide for this node.
  char *MLChildrenChars[4];
  // should only be set on the leaves at first.
  char dna;
} felsensteinData;

char indexToChar(int index) {
  switch(index) {
  case 0: return 'A';
  case 1: return 'G';
  case 2: return 'C';
  case 3: return 'T';
  default: throw hal_exception("Invalid index");
  }
}

int charToIndex(char dna)
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
  default:
    // FIXME: Will die on N's -- need to give some special return value and
    // ignore the char.
    throw hal_exception("Unsupported character " + string(1, dna) + 
                        " found in sequence");
  }
}

// find root genome of tree for this position in the genome
pair<const Genome *, hal_index_t> findRoot(const Genome *genome,
                                           hal_index_t pos) {
  if (genome->getParent() == NULL) {
    return make_pair(genome, pos);
  }
  TopSegmentIteratorConstPtr topIt = genome->getTopSegmentIterator();
  topIt->toSite(pos, false);
  if (!topIt->hasParent()) {
//    cout << "Root genome for pos: " << pos << " " << genome->getName() << endl;
    return make_pair(genome, pos);
  } else {
    BottomSegmentIteratorConstPtr botIt = genome->getParent()->getBottomSegmentIterator();
    botIt->toParent(topIt);
    hal_index_t parentPos = botIt->getStartPosition();
    hal_index_t offset = abs(pos - topIt->getStartPosition());
    parentPos = topIt->getParentReversed() ? parentPos - offset : parentPos + offset;
//    cout << "Going to parent " << genome->getParent()->getName() << ": start pos " << pos << " parent pos: " << parentPos << endl;
    return findRoot(genome->getParent(), parentPos);
  }
}

// Build site-specific tree below this genome.
void buildTree(AlignmentConstPtr alignment, const Genome *genome,
               hal_index_t pos, stTree *tree)
{
  stTree_setLabel(tree, genome->getName().c_str());
  felsensteinData *data = (felsensteinData *) malloc(sizeof(felsensteinData));
  memset(data, 0, sizeof(felsensteinData));
  data->pos = pos;
  stTree_setClientData(tree, data);
  if (genome->getNumChildren() == 0) {
    DNAIteratorConstPtr dnaIt = genome->getDNAIterator(pos);
    data->dna = dnaIt->getChar();
    return;
  }
  for (int i = 0; i < 4; i++) {
    data->MLChildrenChars[i] = (char *) malloc(sizeof(char) * genome->getNumChildren());
    memset(data->MLChildrenChars[i], 0, sizeof(char) * genome->getNumChildren());
  }
  BottomSegmentIteratorConstPtr botIt = genome->getBottomSegmentIterator();
  botIt->toSite(pos, false);
  // Just for debug.
  bool insertion = true;
  for (hal_size_t i = 0; i < botIt->getNumChildren(); i++) {
    hal_index_t childIndex = botIt->getChildIndex(i);
    if (childIndex != NULL_INDEX) {
      insertion = false;
      const Genome *childGenome = genome->getChild(i);
      TopSegmentIteratorConstPtr topIt = childGenome->getTopSegmentIterator();
      topIt->toChild(botIt, i);
      if (topIt->getNextParalogyIndex() != NULL_INDEX) {
        // FIXME: can run into problems comparing these iterators if
        // there is an orientation change between paralogous segments
        // (can possibly create more duplications than really exist)
        // Won't happen with the way the iterator comparison is implemented now
        // though.
        TopSegmentIteratorConstPtr original = topIt->copy();
        for (topIt->toNextParalogy(); !topIt->equals(original); topIt->toNextParalogy()) {
          hal_index_t childPos = topIt->getStartPosition();
          hal_index_t offset = abs(pos - botIt->getStartPosition());
          childPos = botIt->getChildReversed(i) ? childPos - offset : childPos + offset;
          stTree *childNode = stTree_construct();
          stTree_setParent(childNode, tree);
          double branchLength = alignment->getBranchLength(genome->getName(),
                                                           childGenome->getName());
          // FIXME: make option for branch length of duplications
          // (e.g. as fraction of species-tree branch length)
          stTree_setBranchLength(childNode, branchLength);
          buildTree(alignment, childGenome, childPos, childNode);
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
      buildTree(alignment, childGenome, childPos, childNode);
    }
  }
  // Assume that ancestors cannot have insertions. Otherwise will have to think
  // more carefully about how leaves in the stTree are treated.
  assert(insertion == false);
}

double jukesCantor(char childDNA, char parentDNA, double branchLength,
                   double alpha)
{
  if (childDNA == parentDNA) {
    return 0.25*(1 + 3*exp(-4*alpha*branchLength));
  } else {
    return 0.25*(1 - exp(-4*alpha*branchLength));
  }
}

void doFelsenstein(stTree *node, double alpha)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(node);
  if (data->done) {
    return;
  }
  if (stTree_getChildNumber(node) == 0) {
    for (int dna = 0; dna < 4; ++dna) {
      if (dna == charToIndex(data->dna)) {
        cout << stTree_getLabel(node) << ": " << data->dna << endl;
        data->pLeaves[dna] = 1.0;
      } else {
        data->pLeaves[dna] = 0.0;
      }
    }
  } else {
    for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
      doFelsenstein(stTree_getChild(node, childIdx), alpha);
    }


    for (int dna = 0; dna < 4; dna++) {
      double prob = 0.0;
      // Need all possiblities for child nodes, so iterate over the
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
          probAssignment *= jukesCantor(indexToChar(c % 4), indexToChar(dna), stTree_getBranchLength(childNode), alpha);
          probAssignment *= childData->pLeaves[c % 4];
          c /= 4;
        }
//        cout << "probability of assignment " << cartesian << ": " << probAssignment << endl;
        prob += probAssignment;
      }
      data->pLeaves[dna] = prob;
//      cout << stTree_getLabel(node) << ": " << dna << " " << prob << endl;
      // Leave pointers to the letters in each child that are most
      // probable given each letter in this node (could be combined with
      // the above somehow)
      for (int64_t childIdx = 0; childIdx < stTree_getChildNumber(node); childIdx++) {
        stTree *childNode = stTree_getChild(node, childIdx);
        felsensteinData *childData = (felsensteinData *) stTree_getClientData(childNode);
        double maxProb = 0.0;
        int maxDna = -1;
        for (int childDna = 0; childDna < 4; childDna++) {
          double probAssignment = childData->pLeaves[childDna] * jukesCantor(indexToChar(childDna), indexToChar(dna), stTree_getBranchLength(childNode), alpha);
          if (probAssignment > maxProb) {
            maxProb = probAssignment;
            maxDna = childDna;
          }
          if (maxDna == -1) {
            // 0 probability for this assignment, this can happen with
            // leaf nodes since we are sure of their sequence.
            assert(stTree_getChildNumber(childNode) == 0);
            data->MLChildrenChars[dna][childIdx] = '\0';
          } else {
            data->MLChildrenChars[dna][childIdx] = indexToChar(maxDna);
          }
        }
      }
    }
  }
  data->done = true;
}

// Assign nucleotides to each node in the tree.
void walkFelsenstein(stTree *tree, char assignment)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  data->dna = assignment;
  cout << stTree_getLabel(tree) << ": " << assignment << endl;
  for (int i = 0; i < stTree_getChildNumber(tree); i++) {
    stTree *childNode = stTree_getChild(tree, i);
    if (stTree_getChildNumber(childNode) != 0) {
      char childDNA = data->MLChildrenChars[charToIndex(assignment)][i];
      walkFelsenstein(childNode, childDNA);
    }
  }
}

void writeNucleotide(stTree *tree, Genome *targetGenome)
{
  if (targetGenome->getName() == stTree_getLabel(tree)) {
    felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
    DNAIteratorPtr dnaIt = targetGenome->getDNAIterator(data->pos);
    if (data->dna != toupper(dnaIt->getChar())) {
      cout << "different from original char, which is: " << string(1, toupper(dnaIt->getChar())) << " (new char is " << string(1, data->dna) << ")" << endl;
    }
  } else {
    if (stTree_getChildNumber(tree) == 0) {
      return;
    }
    for (int i = 0; i < stTree_getChildNumber(tree); i++) {
      stTree *childNode = stTree_getChild(tree, i);
      writeNucleotide(childNode, targetGenome);
    }
  }
}

void freeClientData(stTree *tree)
{
  felsensteinData *data = (felsensteinData *) stTree_getClientData(tree);
  if (stTree_getChildNumber(tree)) {
    for (int i = 0; i < 4; i++) {
      free(data->MLChildrenChars[i]);
    }
    for (int i = 0; i < stTree_getChildNumber(tree); i++) {
      freeClientData(stTree_getChild(tree, i));
    }
  }
  free(data);
}

int main(int argc, char *argv[])
{
  string halPath, genomeName;
  double alpha;
  CLParserPtr optParser = initParser();
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    genomeName = optParser->getArgument<string>("genome");
    alpha = optParser->getArgument<double>("substRate");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
    throw hal_exception("Genome " + genomeName + " not found in alignment.");
  }
  if (genome->getNumChildren() == 0) {
    throw hal_exception("Genome " + genomeName + " is a leaf genome.");
  }
  for (hal_index_t pos = 0; pos < genome->getSequenceLength(); pos++) {
    stTree *tree = stTree_construct();
    cout << pos << endl;
    // Find root of tree
    pair<const Genome *, hal_index_t> rootInfo = findRoot(genome, pos);
    const Genome *root = rootInfo.first;
    hal_index_t rootPos = rootInfo.second;
    buildTree(alignment, root, rootPos, tree);
    doFelsenstein(tree, alpha);
    // Find assignment for root node that maximizes P(leaves)
    felsensteinData *rootData = (felsensteinData *) stTree_getClientData(tree);
    double maxProb = 0.0;
    int maxDna = -1;
    for (int dna = 0; dna < 4; dna++) {
      if (rootData->pLeaves[dna] > maxProb) {
        maxDna = dna;
        maxProb = rootData->pLeaves[dna];
      }
    }
    walkFelsenstein(tree, indexToChar(maxDna));
    writeNucleotide(tree, genome);
    freeClientData(tree);
    stTree_destruct(tree);
   }
  return 0;
}
