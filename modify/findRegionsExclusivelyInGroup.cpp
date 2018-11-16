#include "hal.h"
#include "halBedLine.h"

using namespace std;
using namespace hal;

static void initParser(CLParser& optionsParser) {
  optionsParser.addArgument("halFile", "hal tree");
  optionsParser.addArgument("referenceGenome", "genome to create the BED file "
                             "for");
  optionsParser.addArgument("ingroupGenomes", "list of 'ingroup' genomes (comma-separated)");
  optionsParser.addOption("minIngroupGenomes", "minimum number of ingroup genomes that a region must appear in (default: all)", -1);
  optionsParser.addOption("maxOutgroupGenomes", "maximum number of outgroup genomes that a region is allowed to be in (default: 0)", 0);
}

int main(int argc, char *argv[])
{
  string halPath, referenceGenomeName, ingroupGenomesUnsplit;
  int minIngroupGenomes, maxOutgroupGenomes;
  CLParser optionsParser;
  initParser(optionsParser);
  try {
    optionsParser.parseOptions(argc, argv);
    halPath = optionsParser.getArgument<string>("halFile");
    referenceGenomeName = optionsParser.getArgument<string>("referenceGenome");
    ingroupGenomesUnsplit = optionsParser.getArgument<string>("ingroupGenomes");
    minIngroupGenomes = optionsParser.getOption<int>("minIngroupGenomes");
    maxOutgroupGenomes = optionsParser.getOption<int>("maxOutgroupGenomes");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optionsParser.printUsage(cerr);
    return 1;
  }
  AlignmentConstPtr alignment(openHalAlignment(halPath, &optionsParser));
  vector<string> ingroupGenomeNames;
  ingroupGenomeNames = chopString(ingroupGenomesUnsplit, ",");
  set <const Genome *>ingroupGenomes;
  for (size_t i = 0; i < ingroupGenomeNames.size(); i++) {
    ingroupGenomes.insert(alignment->openGenome(ingroupGenomeNames[i]));
  }
  // a bit hacky -- -1 is set to the default because it doesn't make sense
  if (minIngroupGenomes == -1) {
    minIngroupGenomes = ingroupGenomes.size();
  }
  const Genome *referenceGenome = alignment->openGenome(referenceGenomeName);
  // FIXME tmp
  ostream &os = cout;
  SequenceIteratorPtr seqIt = referenceGenome->getSequenceIterator();
  for (; not seqIt->atEnd(); seqIt->toNext()) {
    ColumnIteratorPtr colIt = seqIt->getSequence()->getColumnIterator(NULL, 0, 0, NULL_INDEX, false, true);
    bool inRegion = false;
    BedLine curBedLine;
    while (1) {
      bool wasInRegion = inRegion;
      inRegion = true;
      const ColumnIterator::ColumnMap *cols = colIt->getColumnMap();
      ColumnIterator::ColumnMap::const_iterator colMapIt;
      int ingroupCount = 0, outgroupCount = 0;
      // For complicated reasons we can't use the column iterator
      // noDupes setting to ignore duplicates. (It really avoids
      // following paralogies, which can end up causing us to call
      // false positives.) So it's done here instead.
      set <const Genome *> seenGenomes;
      for (colMapIt = cols->begin(); colMapIt != cols->end(); colMapIt++) {
        if (colMapIt->second->empty()) {
          // The column map can contain empty entries.
          continue;
        }
        const Genome *colGenome = colMapIt->first->getGenome();
        if (ingroupGenomes.count(colGenome) && !seenGenomes.count(colGenome)) {
          ingroupCount += 1;
          seenGenomes.insert(colGenome);
        } else if (!seenGenomes.count(colGenome)) {
          outgroupCount += 1;
          seenGenomes.insert(colGenome);
        }
      }
      if (ingroupCount >= minIngroupGenomes && outgroupCount <= maxOutgroupGenomes) {
        inRegion = true;
      } else {
        if (wasInRegion) {
          const hal_index_t pos = colIt->getReferenceSequencePosition();
          curBedLine._end = pos;
          curBedLine.write(os);
        }
        inRegion = false;
      }
      if (inRegion && !wasInRegion) {
        // start of single-copy region, output start of bed entry
        const Sequence * seq = colIt->getReferenceSequence();
        const hal_index_t pos = colIt->getReferenceSequencePosition();
        curBedLine._chrName = seq->getName();
        curBedLine._start = pos;
      }
      if (colIt->lastColumn()) {
        // Have to break here instead of at the beginning of the loop to
        // avoid missing the last column.
        // UPDATE TO LATEST CODE FROM ABOVE (OR STOP BEING LAZY)
        if (inRegion) {
          // current single-copy region has ended, finish bed entry
          const hal_index_t pos = colIt->getReferenceSequencePosition();
          curBedLine._end = pos;
          curBedLine.write(os);
        }
        break;
      }
      colIt->toRight();
    }
  }
  return 0;
}
