#include "hal.h"
#include "ancestorsML.h"
#include "ancestorsMLBed.h"

using namespace std;
using namespace hal;

static CLParserPtr initParser()
{
  CLParserPtr optionsParser = hdf5CLParserInstance(true);
  optionsParser->addArgument("halFile", "hal tree");
  optionsParser->addArgument("genome", "(ancestor) genome to modify");
  optionsParser->addArgument("model", "phyloP model file");
  optionsParser->addOption("startPos", "start position", 0);
  optionsParser->addOption("endPos", "end position", -1);
  optionsParser->addOption("sequence", "Sequence name. IMPORTANT: if "
                           "sequence name is not provided but startPos or "
                           "endPos are, they will be assumed to be in "
                           "genome coordinates", "");
  optionsParser->addOption("thresholdN", "threshold below which an N is output", 0.9);
  optionsParser->addOption("bed", "bed file to scan", "");
  optionsParser->addOptionFlag("outputPosts", "output posterior"
                               " probabilities for reference in wig"
                               " format", false);
  optionsParser->addOptionFlag("printWrites", "print base changes", false);
  return optionsParser;
}


int main(int argc, char *argv[])
{
  string halPath, genomeName, modPath, sequenceName, bedPath;
  CLParserPtr optParser = initParser();
  bool printWrites = false, outputPosts = false;
  hal_index_t startPos = 0;
  hal_index_t endPos = -1;
  double threshold = 0.0;
  try {
    optParser->parseOptions(argc, argv);
    halPath = optParser->getArgument<string>("halFile");
    genomeName = optParser->getArgument<string>("genome");
    modPath = optParser->getArgument<string>("model");
    startPos = optParser->getOption<hal_index_t>("startPos");
    endPos = optParser->getOption<hal_index_t>("endPos");
    threshold = optParser->getOption<double>("thresholdN");
    sequenceName = optParser->getOption<string>("sequence");
    bedPath = optParser->getOption<string>("bed");
    outputPosts = optParser->getFlag("outputPosts");
    printWrites = optParser->getFlag("printWrites");
  } catch (exception &e) {
    optParser->printUsage(cerr);
    return 1;
  }

  // Load phast model.
  FILE *infile = phast_fopen(modPath.c_str(), "r");
  TreeModel *mod = tm_new_from_file(infile, TRUE);
  phast_fclose(infile);
  tm_set_subst_matrices(mod);
  // Map names to phast model IDs.
  map<string, int> nameToId;
  List *phastList = tr_postorder(mod->tree);
  for (int i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = (TreeNode*) lst_get_ptr(phastList, i);
    nameToId[n->name] = n->id;
  }
  lst_free(phastList);

  AlignmentConstPtr alignment = openHalAlignmentReadOnly(halPath, optParser);
  const Genome *genome = alignment->openGenome(genomeName);
  if (genome == NULL) {
    throw hal_exception("Genome " + genomeName + " not found in alignment.");
  }
  if (genome->getNumChildren() == 0) {
    throw hal_exception("Genome " + genomeName + " is a leaf genome.");
  }
  
  if (bedPath != "") {
    AncestorsMLBed bedScanner(mod, alignment, genome, nameToId, threshold,
                              printWrites, outputPosts);
    bedScanner.scan(bedPath, -1);
    return 0;
  }

  if (sequenceName != "") {
    const Sequence *sequence = genome->getSequence(sequenceName);
    if (sequence == NULL) {
      throw hal_exception("Sequence name not found!");
    }
    startPos += sequence->getStartPosition();
    if (endPos == -1) {
      endPos = sequence->getEndPosition();
    } else {
      endPos += sequence->getStartPosition();
      if (endPos > sequence->getEndPosition()) {
        endPos = sequence->getEndPosition();
      }
    }
  }

  if (endPos == -1 || endPos > genome->getSequenceLength()) {
    endPos = genome->getSequenceLength();
  }
  reEstimate(mod, alignment, genome, startPos, endPos, nameToId, threshold, printWrites, outputPosts);
  alignment->close();
  return 0;
}
