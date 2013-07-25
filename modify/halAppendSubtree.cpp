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
  optionsParser->addArgument("bridgeFile", "alignment containing parent,"
                             " subtree root, and its future siblings, if any");
  optionsParser->addArgument("parentName", "node to be added to");
  optionsParser->addArgument("rootName", "name of subtree root");
  optionsParser->addArgument("branchLength", "branch length");
  optionsParser->addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                               " update", false);
  return optionsParser;
}

void addSubtree(AlignmentPtr mainAlignment, AlignmentConstPtr appendAlignment, 
                string currNode)
{
  Genome *outGenome = mainAlignment->openGenome(currNode);
  const Genome *inGenome = appendAlignment->openGenome(currNode);
  inGenome->copy(outGenome);
  vector<string> children = appendAlignment->getChildNames(currNode);
  for (size_t i = 0; i < children.size(); i++) {
    mainAlignment->addLeafGenome(children[i], currNode,
                                 appendAlignment->getBranchLength(currNode,
                                                                   children[i]));
    addSubtree(mainAlignment, appendAlignment, children[i]);
  }
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string mainPath, appendPath, bridgePath, parentName, rootName;
  double branchLength;
  bool noMarkAncestors;
  try {
    optParser->parseOptions(argc, argv);
    mainPath = optParser->getArgument<string>("mainFile");
    appendPath = optParser->getArgument<string>("appendFile");
    bridgePath = optParser->getArgument<string>("bridgeFile");
    parentName = optParser->getArgument<string>("parentName");
    rootName = optParser->getArgument<string>("rootName");
    branchLength = optParser->getArgument<double>("branchLength");
    noMarkAncestors = optParser->getFlag("noMarkAncestors");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment = openHalAlignment(mainPath, optParser);
  AlignmentConstPtr appendAlignment = openHalAlignment(appendPath, optParser);
  AlignmentConstPtr bridgeAlignment = openHalAlignment(bridgePath,
                                                         optParser);

  mainAlignment->addLeafGenome(rootName, parentName, branchLength);
  addSubtree(mainAlignment, appendAlignment, rootName);

  // Need proper bottom segments for parent genome
  Genome *mainParentGenome = mainAlignment->openGenome(parentName);
  const Genome *bridgeParentGenome = bridgeAlignment->openGenome(parentName);
  bridgeParentGenome->copyBottomDimensions(mainParentGenome);
  bridgeParentGenome->copyBottomSegments(mainParentGenome);

  // And top segments for its children
  vector<string> children = bridgeAlignment->getChildNames(parentName);
  for (size_t i = 0; i < children.size(); i++) {
    Genome *mainChildGenome = mainAlignment->openGenome(children[i]);
    const Genome *bridgeChildGenome = bridgeAlignment->openGenome(children[i]);
    bridgeChildGenome->copyTopDimensions(mainChildGenome);
    bridgeChildGenome->copyTopSegments(mainChildGenome);
  }
  if (!noMarkAncestors) {
    markAncestorsForUpdate(mainAlignment, rootName);
  }
  mainAlignment->close();
  appendAlignment->close();
  bridgeAlignment->close();
  return 0;
}
