#include "hal.h"
#include "copyGenomes.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addArgument("botSegmentsFile", "tree containing insert, its "
                             "proper bottom segments, and the new leaf genome");

  // FIXME: now bot, top segments are bad names for the file since
  // there is more copied than just the intermediate top/bot segments:
  // the parent bot segs, child top segs, leaf genome, and
  // intermediate top/bot segs are all copied
  optionsParser->addArgument("topSegmentsFile", "tree containing insert, its "
                             "parent, and its proper top segments");
  optionsParser->addArgument("parentName", "insert's future parent");
  optionsParser->addArgument("insertName", "insert name");
  optionsParser->addArgument("childName", "insert's future child");
  optionsParser->addArgument("leafName", "name of new leaf genome");
  optionsParser->addArgument("upperBranchLength", "length of branch from parent"
                             " to insert");
  optionsParser->addArgument("leafBranchLength", "leaf branch length");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, botSegmentsPath, topSegmentsPath, parentName, insertName,
    childName, leafName;
  double upperBranchLength, leafBranchLength;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    botSegmentsPath = optParser->getArgument<string>("botSegmentsFile");
    topSegmentsPath = optParser->getArgument<string>("topSegmentsFile");
    parentName = optParser->getArgument<string>("parentName");
    insertName = optParser->getArgument<string>("insertName");
    childName = optParser->getArgument<string>("childName");
    leafName = optParser->getArgument<string>("leafName");
    upperBranchLength = optParser->getArgument<double>("upperBranchLength");
    leafBranchLength = optParser->getArgument<double>("leafBranchLength");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment = openHalAlignment(inPath, optParser);
  AlignmentConstPtr botAlignment = openHalAlignment(botSegmentsPath, optParser);
  AlignmentConstPtr topAlignment = openHalAlignment(topSegmentsPath, optParser);
  mainAlignment->insertGenome(insertName, parentName, childName,
                              upperBranchLength);
  mainAlignment->addLeafGenome(leafName, insertName, leafBranchLength);
  // Insert the new intermediate node.
  Genome *insertGenome = mainAlignment->openGenome(insertName);
  const Genome *topInsertGenome = topAlignment->openGenome(insertName);
  const Genome *botInsertGenome = botAlignment->openGenome(insertName);
  copyAllDimensions(topAlignment, topInsertGenome, insertGenome);
  copyTopDimensions(topInsertGenome, insertGenome);
  copyBotDimensions(botInsertGenome, insertGenome);
  copyGenomeWithoutBotSegments(topInsertGenome, insertGenome);
  copyBotSegments(botInsertGenome, insertGenome);
  fixParseInfo(insertGenome);

  // Copy the bottom segments for the parent genome from the top alignment.
  Genome *parentGenome = mainAlignment->openGenome(parentName);
  const Genome *botParentGenome = topAlignment->openGenome(parentName);
  copyBotDimensions(botParentGenome, parentGenome);
  copyBotSegments(botParentGenome, parentGenome);
  fixParseInfo(parentGenome);

  // Fix the parent's other children as well.
  vector<string> allChildren = mainAlignment->getChildNames(parentName);
  for (size_t i = 0; i < allChildren.size(); i++)
  {
    if (allChildren[i] != insertName)
    {
      Genome *outGenome = mainAlignment->openGenome(allChildren[i]);
      const Genome *topSegmentsGenome = topAlignment->openGenome(allChildren[i]);
      copyTopDimensions(topSegmentsGenome, outGenome);
      copyGenomeWithoutBotSegments(topSegmentsGenome, outGenome);
      fixParseInfo(outGenome);
            
    }
  }

  // Copy the top segments for the child genome from the bottom alignment.
  Genome *childGenome = mainAlignment->openGenome(childName);
  const Genome *topChildGenome = botAlignment->openGenome(childName);
  copyTopDimensions(topChildGenome, childGenome);
  copyGenomeWithoutBotSegments(topChildGenome, childGenome);
  fixParseInfo(childGenome);

  // Copy the entire genome for the leaf from the bottom alignment.
  Genome *outLeafGenome = mainAlignment->openGenome(leafName);
  const Genome *inLeafGenome = botAlignment->openGenome(leafName);
  copyAllDimensions(botAlignment, inLeafGenome, outLeafGenome);
  copyGenomeWithoutBotSegments(inLeafGenome, outLeafGenome);

  mainAlignment->close();
  botAlignment->close();
  topAlignment->close();
}
