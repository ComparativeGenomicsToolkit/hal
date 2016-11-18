#include <fstream>
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

}
