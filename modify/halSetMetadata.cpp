// Simple utility to set metadata for a hal genome or alignment
#include "hal.h"

using namespace hal;
using namespace std;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->setDescription("Set metadata for an alignment or genome");
  optionsParser->addArgument("halFile", "hal file to modify");
  optionsParser->addArgument("key", "metadata key");
  optionsParser->addArgument("value", "metadata value");
  optionsParser->addOption("genome", "genome to set metadata for instead of "
                           "setting it for the entire alignment", "");
  return optionsParser;
}

int main(int argc, char *argv[]) {
  string halPath, key, value, genomeName;
  CLParserPtr optParser = initParser();
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    key = optParser->getArgument<string>("key");
    value = optParser->getArgument<string>("value");
    genomeName = optParser->getOption<string>("genome");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }

  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  if (genomeName == "") {
    // No genome to set metadata for, so set the alignment-wide metadata.
    MetaData *metadata = alignment->getMetaData();
    metadata->set(key, value);
  } else {
    Genome *genome = alignment->openGenome(genomeName);
    if (genome == NULL) {
      throw hal_exception("No genome named " + genomeName + " in alignment");
    }
    MetaData *metadata = genome->getMetaData();
    metadata->set(key, value);
  }
  alignment->close();
  return 0;
}
