#include "hal.h"
#include "renameFile.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->setDescription("Rename genomes in a HAL file in-place.");
  optionsParser->addArgument("halFile", "hal file");
  optionsParser->addArgument("renameFile",
                             "Tab-separated file. First column: existing genome"
                             " name, second column: new genome name. Any "
                             "genomes not provided will stay the same.");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string halPath, renamePath;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    renamePath = optParser->getArgument<string>("renameFile");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }

  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  map<string, string> renameMap = ingestRenameFile(renamePath);

  // Check that the alignment has all the old genome names, and none
  // of the new genome names.
  for (map<string, string>::iterator it = renameMap.begin();
       it != renameMap.end(); it++) {
    Genome *genome = alignment->openGenome(it->first);
    if (genome == NULL) {
      throw hal_exception("Genome " + it->first + " not found in alignment");
    }

    genome = alignment->openGenome(it->second);
    if (genome != NULL) {
      throw hal_exception("Attempting to rename " + it->first +
                          " to " + it->second + " failed: " + it->second +
                          " is already in the alignment! Name it to something"
                          " temporary first");
    }
  }

  // Do the actual renaming now that we are relatively sure nothing
  // will go wrong.
  for (map<string, string>::iterator it = renameMap.begin();
       it != renameMap.end(); it++) {
    cout << "Renaming " << it->first << " to " << it->second << endl;
    Genome *genome = alignment->openGenome(it->first);
    genome->rename(it->second);
  }

  alignment->close();
}
