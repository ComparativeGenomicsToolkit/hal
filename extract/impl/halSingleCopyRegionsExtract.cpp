// Creates BED of regions in the specified genome that are in a column
// with no duplicates in any of the target genomes.
#include "hal.h"
#include "halBedLine.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal tree");
  optionsParser->addArgument("referenceGenome", "genome to create the BED file "
                             "for");
  optionsParser->addOption("targetGenomes", "genomes to check for homologous "
                           "duplicated sites (comma-separated, default=leaves)",
                           "");
  optionsParser->addOption("seqStartNum", "sequence number to start with", 0);
  optionsParser->addOption("seqEndNum", "sequence number to end with", -1);
  optionsParser->addOptionFlag("requireAllTargets", "require the regions to be present in all target genomes", false);
  return optionsParser;
}

int main(int argc, char *argv[])
{
  string halPath, referenceGenomeName, targetGenomesUnsplit;
  hal_index_t seqStartPos = 0;
  hal_index_t seqEndPos = -1;
  CLParserPtr optParser = initParser();
  bool requireAllTargets = false;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    referenceGenomeName = optParser->getArgument<string>("referenceGenome");
    targetGenomesUnsplit = optParser->getOption<string>("targetGenomes");
    seqStartPos = optParser->getOption<hal_index_t>("seqStartNum");
    seqEndPos = optParser->getOption<hal_index_t>("seqEndNum");
    requireAllTargets = optParser->getFlag("requireAllTargets");
  } catch (exception &e) {
    cerr << e.what() << endl;
    optParser->printUsage(cerr);
    return 1;
  }
  AlignmentConstPtr alignment = openHalAlignment(halPath, optParser);
  vector<string> targetGenomeNames;
  if (targetGenomesUnsplit != "") {
    // Target genomes provided
    targetGenomeNames = chopString(targetGenomesUnsplit, ",");
  } else {
    // Use all leaf genomes as targets
    targetGenomeNames = alignment->getLeafNamesBelow(alignment->getRootName());
  }
  set<const Genome *> targetGenomes;
  for (size_t i = 0; i < targetGenomeNames.size(); i++) {
    targetGenomes.insert(alignment->openGenome(targetGenomeNames[i]));
  }
  const Genome *referenceGenome = alignment->openGenome(referenceGenomeName);
  ColumnIteratorConstPtr colIt = referenceGenome->getColumnIterator();
  BedLine curBedLine;

  ostream &os = cout;
  SequenceIteratorConstPtr seqIt = referenceGenome->getSequenceIterator(seqStartPos);
  SequenceIteratorConstPtr seqItEnd = (seqEndPos == -1) ? referenceGenome->getSequenceEndIterator() : referenceGenome->getSequenceIterator(seqEndPos + 1);
  for (; seqIt != seqItEnd; seqIt->toNext()) {
    ColumnIteratorConstPtr colIt = seqIt->getSequence()->getColumnIterator();
    bool inRegion = false;
    BedLine curBedLine;
    while (1) {
      bool wasInRegion = inRegion;
      inRegion = true;
      const ColumnIterator::ColumnMap *cols = colIt->getColumnMap();
      ColumnIterator::ColumnMap::const_iterator colMapIt;
      set <const Genome *> seenGenomes;
      unsigned int targetGenomeCount = 0;
      for (colMapIt = cols->begin(); colMapIt != cols->end(); colMapIt++) {
        if (colMapIt->second->empty()) {
          // The column map can contain empty entries.
          continue;
        }
        const Genome *colGenome = colMapIt->first->getGenome();
        if (targetGenomes.count(colGenome)) {
          targetGenomeCount++;
          if (seenGenomes.count(colGenome) || colMapIt->second->size() > 1) {
            // Duplication -- not single-copy among targets.
            inRegion = false;
          }
          seenGenomes.insert(colGenome);
        }    
      }
      if (requireAllTargets && (targetGenomeCount < targetGenomes.size())) {
        // not in every genome
        inRegion = false;
      }
      if (!inRegion && wasInRegion) {
          const hal_index_t pos = colIt->getReferenceSequencePosition();
          curBedLine._end = pos;
          curBedLine.write(os);
      }
      if (inRegion && !wasInRegion) {
        // start of single-copy region, output start of bed entry
        const Sequence * seq = seqIt->getSequence();
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
      if (colIt->getReferenceSequencePosition() % 10000 == 0) {
        colIt->defragment();
      }
      // colIt->toRight();
      colIt->toSite(colIt->getReferenceSequencePosition() + seqIt->getSequence()->getStartPosition() + 1, seqIt->getSequence()->getEndPosition(), true);
    }
  }
}
