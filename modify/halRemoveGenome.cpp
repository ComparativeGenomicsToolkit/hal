#include "hal.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("inFile", "existing tree");
  optionsParser->addArgument("deleteNode", "(leaf) genome to delete");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  optParser->parseOptions(argc, argv);
  string inPath = optParser->getArgument<string>("inFile");
  string deleteNode = optParser->getArgument<string>("deleteNode");
  AlignmentPtr alignment = openHalAlignment(inPath, optParser);
  alignment->removeGenome(deleteNode);
  return 0;
}
