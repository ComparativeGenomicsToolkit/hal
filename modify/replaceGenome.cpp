// Need to replace a genome with its top & bottom segments. (2 files needed)
#include "hal.h"
#include "copyGenomes.h"

using namespace hal;
using namespace std;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addArgument("botAlignmentFile", "hal file containing an "
                             "alignment of the genome and its children");
  optionsParser->addArgument("topAlignmentFile", "hal file containing an "
                             "alignment of the genome, its parent, and "
                             "its siblings");
  optionsParser->addArgument("genomeName", "name of genome to be replaced");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, botAlignmentFile, topAlignmentFile, genomeName;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    botAlignmentFile = optParser->getArgument<string>("botAlignmentFile");
    topAlignmentFile = optParser->getArgument<string>("topAlignmentFile");
    genomeName = optParser->getArgument<string>("genomeName");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr mainAlignment = openHalAlignment(inPath, optParser);
  AlignmentConstPtr topAlignment = openHalAlignment(topAlignmentFile,
                                                        optParser);
  AlignmentConstPtr botAlignment = openHalAlignment(botAlignmentFile,
                                                        optParser);

  // Copy genome & segments for the genome that's being replaced
  Genome *mainReplacedGenome = mainAlignment->openGenome(genomeName);
  const Genome *topReplacedGenome = topAlignment->openGenome(genomeName);
  const Genome *botReplacedGenome = botAlignment->openGenome(genomeName);
  copyAllDimensions(topAlignment, topReplacedGenome, mainReplacedGenome);
  copyTopDimensions(topReplacedGenome, mainReplacedGenome);
  copyBotDimensions(botReplacedGenome, mainReplacedGenome);
  copyGenomeWithoutBotSegments(topReplacedGenome, mainReplacedGenome);
  copyBotSegments(botReplacedGenome, mainReplacedGenome);
  fixParseInfo(mainReplacedGenome);

  // Copy bot segments for the parent and top segments for the
  // siblings of the genome that's being replaced
  Genome *mainParent = mainReplacedGenome->getParent();
  const Genome *topParent = topReplacedGenome->getParent();
  copyBotDimensions(topParent, mainParent);
  copyBotSegments(topParent, mainParent);
  vector<string> parentChildren = mainAlignment->getChildNames(mainParent->getName());
  for (size_t i = 0; i < parentChildren.size(); i++)
  {
    if (parentChildren[i] != genomeName)
    {
      Genome *mainChild = mainAlignment->openGenome(parentChildren[i]);
      const Genome *topChild  = topAlignment->openGenome(parentChildren[i]);
      copyTopDimensions(topChild, mainChild);
      copyGenomeWithoutBotSegments(topChild, mainChild);
      fixParseInfo(mainChild);
    }
  }
  // Copy top segments for the children
  vector<string> children = mainAlignment->getChildNames(genomeName);  
  for (size_t i = 0; i < children.size(); i++)
  {
    Genome *mainChild = mainAlignment->openGenome(children[i]);
    const Genome *topChild  = topAlignment->openGenome(parentChildren[i]);
    copyTopDimensions(topChild, mainChild);
    copyGenomeWithoutBotSegments(topChild, mainChild);
    fixParseInfo(mainChild);
  }

  mainAlignment->close();
  botAlignment->close();
  topAlignment->close();
}
