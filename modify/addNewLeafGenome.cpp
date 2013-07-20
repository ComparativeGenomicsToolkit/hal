#include "hal.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addArgument("subtreeAlignmentFile", "hal file containing an "
                             "alignment of the parent, its children, and the "
                             "new leaf");
  optionsParser->addArgument("parentName", "parent of new leaf");
  optionsParser->addArgument("leafName", "name of new leaf genome");
  optionsParser->addArgument("leafBranchLength", "leaf branch length");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, subtreeAlignmentFile, parentName, leafName;
  double leafBranchLength;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    subtreeAlignmentFile = optParser->getArgument<string>("subtreeAlignmentFile");
    parentName = optParser->getArgument<string>("parentName");
    leafName = optParser->getArgument<string>("leafName");
    leafBranchLength = optParser->getArgument<double>("leafBranchLength");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }

  AlignmentPtr mainAlignment = openHalAlignment(inPath, optParser);
  AlignmentConstPtr subtreeAlignment = openHalAlignment(subtreeAlignmentFile,
                                                        optParser);
    
  Genome *outLeafGenome = mainAlignment->addLeafGenome(leafName, parentName,
                                                       leafBranchLength);
  const Genome *inLeafGenome = subtreeAlignment->openGenome(leafName);
  inLeafGenome->copy(outLeafGenome);
    
  Genome *outParentGenome = mainAlignment->openGenome(parentName);
  const Genome *inParentGenome = subtreeAlignment->openGenome(parentName);
  inParentGenome->copyBottomDimensions(outParentGenome);
  inParentGenome->copyBottomSegments(outParentGenome);
  outParentGenome->fixParseInfo();

  // Fix the parent's other children as well.
  vector<string> allChildren = mainAlignment->getChildNames(parentName);
  for (size_t i = 0; i < allChildren.size(); i++)
  {
    if (allChildren[i] != leafName)
    {
      Genome *outGenome = mainAlignment->openGenome(allChildren[i]);
      const Genome *inGenome = subtreeAlignment->openGenome(allChildren[i]);
      inGenome->copyTopDimensions(outGenome);
      inGenome->copyTopSegments(outGenome);
      outGenome->fixParseInfo();
    }
  }
}
