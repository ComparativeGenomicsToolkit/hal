#include "hal.h"
#include "markAncestors.h"

using namespace hal;
using namespace std;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addOption("bottomAlignmentFile", "hal file containing an "
                           "alignment of the genome and its children. "
                           "Required for non-leaf genomes.", "\"\"");
  optionsParser->addOption("topAlignmentFile", "hal file containing an "
                           "alignment of the genome, its parent, and "
                           "its siblings. Required if the genome to be "
                           "replaced is not the root.", "\"\"");
  optionsParser->addArgument("genomeName", "name of genome to be replaced");
  optionsParser->addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                               " update", false);
  return optionsParser;
}

void copyFromTopAlignment(AlignmentConstPtr topAlignment,
                          AlignmentPtr mainAlignment, const string &genomeName)
{
  Genome *mainReplacedGenome = mainAlignment->openGenome(genomeName);
  const Genome *topReplacedGenome = topAlignment->openGenome(genomeName);
  topReplacedGenome->copyTopDimensions(mainReplacedGenome);
  topReplacedGenome->copyTopSegments(mainReplacedGenome);
  mainReplacedGenome->fixParseInfo();
  // Copy bot segments for the parent and top segments for the
  // siblings of the genome that's being replaced
  Genome *mainParent = mainReplacedGenome->getParent();
  const Genome *topParent = topReplacedGenome->getParent();
  topParent->copyBottomDimensions(mainParent);
  topParent->copyBottomSegments(mainParent);
  mainParent->fixParseInfo();
  vector<string> siblings = mainAlignment->getChildNames(mainParent->getName());
  for (size_t i = 0; i < siblings.size(); i++)
  {
    if (siblings[i] != genomeName)
    {
      Genome *mainChild = mainAlignment->openGenome(siblings[i]);
      const Genome *topChild  = topAlignment->openGenome(siblings[i]);
      topChild->copyTopDimensions(mainChild);
      topChild->copyTopSegments(mainChild);
      mainChild->fixParseInfo();
    }
  }
}

void copyFromBottomAlignment(AlignmentConstPtr bottomAlignment,
                             AlignmentPtr mainAlignment,
                             const string &genomeName)
{
  // Copy genome & bottom segments for the genome that's being replaced
  Genome *mainReplacedGenome = mainAlignment->openGenome(genomeName);
  const Genome *botReplacedGenome = bottomAlignment->openGenome(genomeName);
  botReplacedGenome->copyDimensions(mainReplacedGenome);
  botReplacedGenome->copySequence(mainReplacedGenome);
  botReplacedGenome->copyBottomDimensions(mainReplacedGenome);
  botReplacedGenome->copyBottomSegments(mainReplacedGenome);
  mainReplacedGenome->fixParseInfo();

  // Copy top segments for the children
  vector<string> children = mainAlignment->getChildNames(genomeName);  
  for (size_t i = 0; i < children.size(); i++)
  {
    Genome *mainChild = mainAlignment->openGenome(children[i]);
    const Genome *botChild  = bottomAlignment->openGenome(children[i]);
    botChild->copyTopDimensions(mainChild);
    botChild->copyTopSegments(mainChild);
    mainChild->fixParseInfo();
  }
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, bottomAlignmentFile, topAlignmentFile, genomeName;
  bool noMarkAncestors;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    bottomAlignmentFile = optParser->getOption<string>("bottomAlignmentFile");
    topAlignmentFile = optParser->getOption<string>("topAlignmentFile");
    genomeName = optParser->getArgument<string>("genomeName");
    noMarkAncestors = optParser->getFlag("noMarkAncestors");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment = openHalAlignment(inPath, optParser);
  AlignmentConstPtr bottomAlignment;
  AlignmentConstPtr topAlignment;
  bool useTopAlignment = mainAlignment->getRootName() != genomeName;
  bool useBottomAlignment = mainAlignment->getChildNames(genomeName).size() != 0;
  Genome *mainReplacedGenome = mainAlignment->openGenome(genomeName);
  if (useTopAlignment) {
    // Not a root genome. Can update using a top alignment.
    if (topAlignmentFile == "\"\"") {
      throw hal_exception("Cannot replace non-root genome without a top "
                          "alignment file.");
    }
    topAlignment = openHalAlignment(topAlignmentFile,
                                    optParser);
    const Genome *topReplacedGenome = topAlignment->openGenome(genomeName);
    topReplacedGenome->copyDimensions(mainReplacedGenome);
    topReplacedGenome->copySequence(mainReplacedGenome);
    
  }
  if (useBottomAlignment) {
    // Not a leaf genome. Can update using a bottom alignment.
    if (bottomAlignmentFile == "\"\"") {
      throw hal_exception("Cannot replace non-leaf genome without a bottom "
                          "alignment file.");
    }
    bottomAlignment = openHalAlignment(bottomAlignmentFile, optParser);
    const Genome *botReplacedGenome = bottomAlignment->openGenome(genomeName);
    botReplacedGenome->copyDimensions(mainReplacedGenome);
    botReplacedGenome->copySequence(mainReplacedGenome);
  }
  if (!useTopAlignment && !useBottomAlignment) {
    throw hal_exception("Root genome is also a leaf genome.");
  }
  if (useBottomAlignment) {
    copyFromBottomAlignment(bottomAlignment, mainAlignment, genomeName);
  }
  if (useTopAlignment) {
    copyFromTopAlignment(topAlignment, mainAlignment, genomeName);
  }

  // Clear update flag if present, since the genome has just been updated.
  MetaData *metaData = mainReplacedGenome->getMetaData();
  if (metaData->has("needsUpdate")) {
    metaData->set("needsUpdate", "false");
  }

  if (!noMarkAncestors) {
    markAncestorsForUpdate(mainAlignment, genomeName);
  }
  if (useTopAlignment) {
    topAlignment->close();
  }
  if (useBottomAlignment) {
    bottomAlignment->close();
  }
  mainAlignment->close();
}
