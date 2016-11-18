#include <fstream>
#include "hal.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal file");
  optionsParser->addArgument("renameFile",
                             "Tab-separated file. First column: existing genome"
                             " name, second column: new genome name");
  return optionsParser;
}

static map<string, string> ingestRenameFile(string tsvPath) {
  map<string, string> ret;
  set<string> values;
  ifstream tsv(tsvPath.c_str());
  string line;
  while (getline(tsv, line)) {
    stringstream lineStream(line);
    string oldName;
    string newName;
    lineStream >> oldName;
    lineStream >> newName;

    if (ret.count(oldName)) {
      throw hal_exception("Old name " + oldName +
                          " is represented twice in rename file");
    }
    if (values.count(oldName)) {
      throw hal_exception("New name " + oldName + " is same as an old name " +
                          oldName + ". Collisions are not allowed, rename to "
                          "something temporary first." );
    }
    if (values.count(newName)) {
      throw hal_exception("New name " + newName +
                          " is represented twice in rename file");
    }

    ret.insert(make_pair(oldName, newName));
    values.insert(newName);
  }
  return ret;
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
