#include "hal.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addArgument("deleteNode", "(leaf) genome to delete");
  optionsParser->addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                               " update", false);
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, deleteNode;
  bool noMarkAncestors;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    deleteNode = optParser->getArgument<string>("deleteNode");
    noMarkAncestors = optParser->getFlag("noMarkAncestors");
  } catch (exception &e) {
    optParser->printUsage(cerr);
  }
  AlignmentPtr alignment = openHalAlignment(inPath, optParser);
  if (!noMarkAncestors) {
    markAncestorsForUpdate(alignment, deleteNode);
  }
  alignment->removeGenome(deleteNode);
  return 0;
}
