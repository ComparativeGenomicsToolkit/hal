#include "hal.h"
#include "sonLibTree.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal file");
  optionsParser->addArgument("newickTree", "newick tree (must be identical,"
                             " except for the branch lengths");
  return optionsParser;
}

void updateBranches(AlignmentPtr alignment, Genome *genome, stTree *newTree)
{
  if (genome->getNumChildren() == 0) {
    return;
  }
  for (hal_size_t i = 0; i < genome->getNumChildren(); i++) {
    Genome *child = genome->getChild(i);
    bool found = false;
    for (int64_t j = 0; j < stTree_getChildNumber(newTree); j++) {
      stTree *newChild = stTree_getChild(newTree, j);
      if (child->getName() == stTree_getLabel(newChild)) {
        found = true;
        alignment->updateBranchLength(genome->getName(), child->getName(),
                                      stTree_getBranchLength(newChild));
        updateBranches(alignment, child, newChild);
        break;
      }
    }
    if (!found) {
      throw hal_exception("Genome " + child->getName() + " not found in proper"
                          " place in replacement newick tree.");
    }
  }
}

int main(int argc, char *argv[])
{
  string halPath, newickTree;
  CLParserPtr optParser = initParser();
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    newickTree = optParser->getArgument<string>("newickTree");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  stTree *newTree = stTree_parseNewickString(newickTree.c_str());
  // recursively update branches
  updateBranches(alignment,
                 alignment->openGenome(alignment->getRootName()), newTree);
}
