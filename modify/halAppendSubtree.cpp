#include "hal.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("mainFile", "destination tree");
  optionsParser->addArgument("appendFile", "alignment containing the tree to be"
                             " appended");
  optionsParser->addArgument("parentName", "node to be added to");
  optionsParser->addArgument("rootName", "name of subtree root");
  optionsParser->addOption("bridgeFile", "alignment containing parent,"
                           " subtree root, and its future siblings, if any "
                           "(required if not merging appended and appendee "
                           "nodes)", "");
  optionsParser->addOption("branchLength", "branch length between appended and "
                           "appendee nodes", 0.0);
  optionsParser->addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                               " update", false);
  optionsParser->addOptionFlag("merge", "merge appended root and node that is appended to",
                               false);
  return optionsParser;
}

void addSubtree(AlignmentPtr mainAlignment, AlignmentConstPtr appendAlignment, 
                string currNode)
{
  Genome *outGenome = mainAlignment->openGenome(currNode);
  const Genome *inGenome = appendAlignment->openGenome(currNode);
  vector<string> children = appendAlignment->getChildNames(currNode);
  for (size_t i = 0; i < children.size(); i++) {
    if (mainAlignment->openGenome(children[i]) != NULL) {
      // Special case for a node that is being merged; we don't want
      // to append children that already exist.
      continue;
    }
    Genome *mainChildGenome = mainAlignment->addLeafGenome(children[i], currNode, appendAlignment->getBranchLength(currNode, children[i]));
    const Genome *appendChildGenome = appendAlignment->openGenome(children[i]);
    appendChildGenome->copy(mainChildGenome);
    addSubtree(mainAlignment, appendAlignment, children[i]);
  }
  inGenome->copyBottomDimensions(outGenome);
  inGenome->copyBottomSegments(outGenome);
  outGenome->fixParseInfo();
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string mainPath, appendPath, bridgePath, parentName, rootName;
  double branchLength;
  bool noMarkAncestors;
  bool merge;
  try {
    optParser->parseOptions(argc, argv);
    mainPath = optParser->getArgument<string>("mainFile");
    appendPath = optParser->getArgument<string>("appendFile");
    bridgePath = optParser->getOption<string>("bridgeFile");
    parentName = optParser->getArgument<string>("parentName");
    rootName = optParser->getArgument<string>("rootName");
    branchLength = optParser->getOption<double>("branchLength");
    noMarkAncestors = optParser->getFlag("noMarkAncestors");
    merge = optParser->getFlag("merge");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment = openHalAlignment(mainPath, optParser);
  AlignmentConstPtr appendAlignment = openHalAlignment(appendPath, optParser);
  AlignmentConstPtr bridgeAlignment;

  if (!merge) {
    if (bridgePath == "") {
      throw hal_exception("need a bridge alignment if not merging nodes");
    }
    bridgeAlignment = openHalAlignment(bridgePath, optParser);
    Genome *mainAppendedRoot = mainAlignment->addLeafGenome(rootName,
                                                            parentName,
                                                            branchLength);
    const Genome *appendAppendedRoot = appendAlignment->openGenome(rootName);
    appendAppendedRoot->copy(mainAppendedRoot);
  } else {
    // the bridge alignment is equivalent to the append alignment in this case
    // (the append alignment will contain at least all the information that
    // the bridge alignment would)
    bridgeAlignment = appendAlignment;
    if (parentName != rootName) {
      throw hal_exception("parent name must be equal to root name if "
                          "--merge option is given");
    }
    assert(branchLength == 0.0);
  }
  addSubtree(mainAlignment, appendAlignment, rootName);

  // Need proper bottom segments for parent genome
  Genome *mainParentGenome = mainAlignment->openGenome(parentName);
  const Genome *bridgeParentGenome = bridgeAlignment->openGenome(parentName);
  bridgeParentGenome->copyBottomDimensions(mainParentGenome);
  bridgeParentGenome->copyBottomSegments(mainParentGenome);
  mainParentGenome->fixParseInfo();

  // And top segments for its children
  vector<string> children = bridgeAlignment->getChildNames(parentName);
  for (size_t i = 0; i < children.size(); i++) {
    Genome *mainChildGenome = mainAlignment->openGenome(children[i]);
    const Genome *bridgeChildGenome = bridgeAlignment->openGenome(children[i]);
    bridgeChildGenome->copyTopDimensions(mainChildGenome);
    bridgeChildGenome->copyTopSegments(mainChildGenome);
    mainChildGenome->fixParseInfo();
  }

  if (!noMarkAncestors) {
    markAncestorsForUpdate(mainAlignment, rootName);
  }
  mainAlignment->close();
  appendAlignment->close();
  if (!merge) {
    bridgeAlignment->close();
  }
  return 0;
}
