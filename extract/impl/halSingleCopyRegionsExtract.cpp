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
  optionsParser->addOption("refSequence", "sequence to traverse", "");
  optionsParser->addOption("start", "start position within the sequence "
                           "(within entire genome if --refSequence is not "
                           "set)", 0);
  optionsParser->addOption("length", "length to traverse (default: until end "
                           "of genome/sequence)", -1);
  optionsParser->addOptionFlag("requireAllTargets", "require the regions to be present in all target genomes", false);
  return optionsParser;
}

int main(int argc, char *argv[])
{
  string halPath, referenceGenomeName, targetGenomesUnsplit, refSequence;
  hal_index_t start, length;
  CLParserPtr optParser = initParser();
  bool requireAllTargets = false;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    referenceGenomeName = optParser->getArgument<string>("referenceGenome");
    targetGenomesUnsplit = optParser->getOption<string>("targetGenomes");
    refSequence = optParser->getOption<string>("refSequence");
    start = optParser->getOption<hal_index_t>("start");
    length = optParser->getOption<hal_index_t>("length");
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
  if (referenceGenome == NULL) {
    throw hal_exception("Genome " + referenceGenomeName + " not present in alignment");
  }

  ostream &os = cout;
  const SegmentedSequence *sequence;
  size_t seqStart, seqEnd;
  if (refSequence != "") {
    sequence = referenceGenome->getSequence(refSequence);
    seqStart = referenceGenome->getSequence(refSequence)->getStartPosition();
    seqEnd = referenceGenome->getSequence(refSequence)->getEndPosition();
  } else {
    sequence = referenceGenome;
    seqStart = 0;
    seqEnd = referenceGenome->getSequenceLength() - 1;
  }

  if (length > (hal_index_t) (seqEnd - seqStart - start + 1)) {
    throw hal_exception("region too long, goes off the end of sequence or genome.");
  }

  ColumnIteratorConstPtr colIt = referenceGenome->getColumnIterator(&targetGenomes,
                                                                    0,
                                                                    start + seqStart,
                                                                    length == -1 ? seqEnd : seqStart + start + length - 1);
  bool inRegion = false;
  BedLine curBedLine;
  const Sequence *prevSequence = NULL;
  hal_index_t prevPos = NULL_INDEX;
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
    bool changedSequence = prevSequence != NULL && prevSequence != colIt->getReferenceSequence();
    if (inRegion && !wasInRegion) {
      // start of single-copy region, output start of bed entry
      const Sequence *seq = colIt->getReferenceSequence();
      const hal_index_t pos = colIt->getReferenceSequencePosition();
      curBedLine._chrName = seq->getName();
      curBedLine._start = pos;
    }
    if (changedSequence || colIt->lastColumn()) {
      if (colIt->lastColumn() && inRegion && length == 1) {
        curBedLine._chrName = colIt->getReferenceSequence()->getName();
        curBedLine._start = colIt->getReferenceSequencePosition();
        curBedLine._end = colIt->getReferenceSequencePosition() + 1;
        curBedLine.write(os);
      } else if (wasInRegion) {
        // current single-copy region has ended, finish bed entry
        curBedLine._end = changedSequence ? prevPos + 1 : colIt->getReferenceSequencePosition() + 1;
        curBedLine.write(os);
        if (inRegion) {
          curBedLine._chrName = colIt->getReferenceSequence()->getName();
          curBedLine._start = 0;
        }
      }
      if (colIt->lastColumn()) {
        // Have to break here instead of at the beginning of the loop to
        // avoid missing the last column.
        break;
      }
    } else if (!inRegion && wasInRegion) {
      const hal_index_t pos = colIt->getReferenceSequencePosition();
      curBedLine._end = pos;
      curBedLine.write(os);
    }
    if (colIt->getReferenceSequencePosition() % 10000 == 0) {
      colIt->defragment();
    }
    prevSequence = colIt->getReferenceSequence();
    prevPos = colIt->getReferenceSequencePosition();
    colIt->toSite(colIt->getReferenceSequencePosition() + colIt->getReferenceSequence()->getStartPosition() + 1, length == -1 ? seqEnd : seqStart + start + length - 1, true);
  }
}
