#include "hal.h"
#include "markAncestors.h"
#include "halCLParser.h"
#include "halAlignmentInstance.h"

using namespace std;
using namespace hal;

static void initParser(CLParser& optionsParser) {
  optionsParser.addArgument("mainFile", "destination tree");
  optionsParser.addArgument("appendFile", "alignment containing the tree to be"
                             " appended");
  optionsParser.addArgument("parentName", "node to be added to");
  optionsParser.addArgument("rootName", "name of subtree root");
  optionsParser.addOption("bridgeFile", "alignment containing parent,"
                           " subtree root, and its future siblings, if any "
                           "(required if not merging appended and appendee "
                           "nodes)", "");
  optionsParser.addOption("branchLength", "branch length between appended and "
                           "appendee nodes", 0.0);
  optionsParser.addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                               " update", false);
  optionsParser.addOptionFlag("merge", "merge appended root and node that is appended to",
                               false);
}

void addSubtree(Alignment* mainAlignment, const Alignment* appendAlignment, 
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
  CLParser optionsParser(WRITE_ACCESS);
  initParser(optionsParser);
  string mainPath, appendPath, bridgePath, parentName, rootName;
  double branchLength;
  bool noMarkAncestors;
  bool merge;
  try {
    optionsParser.parseOptions(argc, argv);
    mainPath = optionsParser.getArgument<string>("mainFile");
    appendPath = optionsParser.getArgument<string>("appendFile");
    bridgePath = optionsParser.getOption<string>("bridgeFile");
    parentName = optionsParser.getArgument<string>("parentName");
    rootName = optionsParser.getArgument<string>("rootName");
    branchLength = optionsParser.getOption<double>("branchLength");
    noMarkAncestors = optionsParser.getFlag("noMarkAncestors");
    merge = optionsParser.getFlag("merge");
  } catch (exception &e) {
    optionsParser.printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment(openHalAlignment(mainPath, &optionsParser));
  AlignmentConstPtr appendAlignment(openHalAlignment(appendPath, &optionsParser));
  AlignmentConstPtr bridgeAlignment;

  if (!merge) {
    if (bridgePath == "") {
      throw hal_exception("need a bridge alignment if not merging nodes");
    }
    bridgeAlignment = AlignmentConstPtr(openHalAlignment(bridgePath, &optionsParser));
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
  addSubtree(mainAlignment.get(), appendAlignment.get(), rootName);

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
    markAncestorsForUpdate(mainAlignment.get(), rootName);
  }
  return 0;
}
