#include "hal.h"
#include "renameFile.h"

using namespace hal;
using namespace std;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->setDescription("Rename the sequences of a genome in-place.");
  optionsParser->addArgument("halFile", "hal file");
  optionsParser->addArgument("genome", "genome to rename the sequences of");
  optionsParser->addArgument("renameFile",
                             "Tab-separated file. First column: existing "
                             "sequence name, second column: new sequence name."
                             " Any sequences not provided will stay the same.");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string halPath, renamePath, genomeName;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    genomeName = optParser->getArgument<string>("genome");
    renamePath = optParser->getArgument<string>("renameFile");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }

  AlignmentPtr alignment = openHalAlignment(halPath, optParser);
  Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
      throw hal_exception("Genome " + genomeName + " not found in alignment");
  }
  map<string, string> renameMap = ingestRenameFile(renamePath);

  for (map<string, string>::iterator it = renameMap.begin();
       it != renameMap.end(); it++) {
    Sequence *sequence = genome->getSequence(it->first);
    if (sequence == NULL) {
      throw hal_exception("Sequence " + it->first + " not found in alignment");
    }

    sequence = genome->getSequence(it->second);
    if (sequence != NULL) {
      throw hal_exception("Attempting to rename sequence " + it->first +
                          " to " + it->second + " failed: " + it->second +
                          " is already in the genome! Name it to something"
                          " temporary first");
    }
  }

  // Do the actual renaming now that we are relatively sure nothing
  // will go wrong.
  for (map<string, string>::iterator it = renameMap.begin();
       it != renameMap.end(); it++) {
    cout << "Renaming " << it->first << " to " << it->second << endl;
    Sequence *sequence = genome->getSequence(it->first);
    sequence->setName(it->second);
  }

  alignment->close();
}
