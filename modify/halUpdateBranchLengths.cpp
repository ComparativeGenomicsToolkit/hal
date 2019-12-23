#include "hal.h"
#include "halAlignmentInstance.h"
#include "sonLibTree.h"

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("halFile", "hal file");
    optionsParser.addArgument("newickTree", "newick tree (must be identical,"
                                            " except for the branch lengths");
}

void updateBranches(Alignment *alignment, Genome *genome, stTree *newTree) {
    if (genome->getNumChildren() == 0) {
        return;
    }
    for (hal_size_t i = 0; i < genome->getNumChildren(); i++) {
        Genome *child = genome->getChild(i);
        bool found = false;
        for (int64_t j = 0; j < stTree_getChildNumber(newTree); j++) {
            stTree *newChild = stTree_getChild(newTree, j);
            if (child->getName() == stTree_getLabel(newChild)) {
                found = true;
                alignment->updateBranchLength(genome->getName(), child->getName(), stTree_getBranchLength(newChild));
                updateBranches(alignment, child, newChild);
                break;
            }
        }
        if (!found) {
            throw hal_exception("Genome " + child->getName() + " not found in proper"
                                                               " place in replacement newick tree.");
        }
    }
}

int main(int argc, char *argv[]) {
    string halPath, newickTree;
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        newickTree = optionsParser.getArgument<string>("newickTree");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        return 1;
    }
    AlignmentPtr alignment(openHalAlignment(halPath, &optionsParser,
                                            WRITE_ACCESS | READ_ACCESS));
    stTree *newTree = stTree_parseNewickString(newickTree.c_str());
    // recursively update branches
    updateBranches(alignment.get(), alignment->openGenome(alignment->getRootName()), newTree);
}
