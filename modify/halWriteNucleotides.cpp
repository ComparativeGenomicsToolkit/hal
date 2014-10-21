// Hacky script to allow writing to a genome's sequence from a TSV (probably
// after an ancestorsML run).
#include <fstream>
#include "hal.h"

using namespace hal;
using namespace std;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->setDescription("Write changes to a hal sequence from a TSV "
                                "containing fields "
                                "genomeName\tpos\toldChar\tnewChar. Note that "
                                "the position is in genome coordinates!");
  optionsParser->addArgument("inFile", "hal file");
  optionsParser->addArgument("tsvFile", "tsv file");
  return optionsParser;
}

int main(int argc, char *argv[])
{
  CLParserPtr optParser = initParser();
  string inPath, tsvFile;
  try {
    optParser->parseOptions(argc, argv);
    inPath = optParser->getArgument<string>("inFile");
    tsvFile = optParser->getArgument<string>("tsvFile");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentPtr alignment = openHalAlignment(inPath, optParser);

  ifstream tsv(tsvFile.c_str());
  string line;
  int64_t lineNum = 0;
  while (getline(tsv, line)) {
    stringstream lineStream(line);
    string genomeName;
    hal_index_t pos;
    char prevChar, newChar;

    lineNum++;
    if (lineNum % 100000 == 0) {
      cout << lineNum << endl;
    }
    lineStream >> genomeName;
    lineStream >> pos;
    lineStream >> prevChar;
    lineStream >> newChar;
    Genome *genome = alignment->openGenome(genomeName);
    DNAIteratorPtr dnaIt = genome->getDNAIterator(pos);
    if (toupper(dnaIt->getChar()) != prevChar) {
      dnaIt->toReverse();
      if (toupper(dnaIt->getChar()) != prevChar) {
        throw hal_exception("previous nucleotide " + string(1, dnaIt->getChar()) + " does not match expected " + string(1, prevChar) + "! Aborting early. Your hal file could be invalid now.");
      }
    }

    dnaIt->setChar(newChar);
  }
  tsv.close();
  alignment->close();
  return 0;
}
