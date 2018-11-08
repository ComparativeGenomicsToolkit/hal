#include "hal.h"
#include "renameFile.h"

using namespace std;
using namespace hal;

static void initParser(CLParser& optionsParser) {
  optionsParser.setDescription("Rename genomes in a HAL file in-place.");
  optionsParser.addArgument("halFile", "hal file");
  optionsParser.addArgument("renameFile",
                             "Tab-separated file. First column: existing genome"
                             " name, second column: new genome name. Any "
                             "genomes not provided will stay the same.");
}

int main(int argc, char *argv[])
{
  CLParser optionsParser(WRITE_ACCESS);
  initParser(optionsParser);
  string halPath, renamePath;
  try {
    optionsParser.parseOptions(argc, argv);
    halPath = optionsParser.getArgument<string>("halFile");
    renamePath = optionsParser.getArgument<string>("renameFile");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optionsParser.printUsage(cerr);
    return 1;
  }

  AlignmentPtr alignment(openHalAlignment(halPath, &optionsParser));
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
