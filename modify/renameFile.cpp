#include <fstream>
#include "sonLib.h"
#include "hal.h"

using namespace std;
using namespace hal;

namespace hal {

map<string, string> ingestRenameFile(string tsvPath) {
  map<string, string> ret;
  set<string> values;
  ifstream tsv(tsvPath.c_str());
  string line;
  while (getline(tsv, line)) {
    stList *tokens = stString_splitByString(line.c_str(), "\t");
    if (stList_length(tokens) != 2) {
      throw hal_exception("Rename file does not have 2 tab-separated fields "
                          "in line: " + line);
    }
    string oldName((char *) stList_get(tokens, 0));
    string newName((char *) stList_get(tokens, 1));
    stList_destruct(tokens);

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

}
