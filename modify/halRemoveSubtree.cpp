#include "hal.h"
#include "halCLParser.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("inFile", "existing tree");
    optionsParser.addArgument("root", "subtree below this node will be deleted (but not the node itself)");
    optionsParser.addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                                                   " update",
                                false);
}

void remove_recursive(AlignmentPtr aln, const string& node) {
    vector<string> children = aln->getChildNames(node);
    for (size_t i = 0; i < children.size(); ++i) {
        remove_recursive(aln, children[i]);
    }
    cerr << "[halRemoveSubtree] removing " << node << endl;
    aln->removeGenome(node);    
}

int main(int argc, char *argv[]) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    string inPath, root;
    bool noMarkAncestors;
    try {
        optionsParser.parseOptions(argc, argv);
        inPath = optionsParser.getArgument<string>("inFile");
        root = optionsParser.getArgument<string>("root");
        noMarkAncestors = optionsParser.getFlag("noMarkAncestors");
    } catch (exception &e) {
        optionsParser.printUsage(cerr);
        return 1;
    }
    AlignmentPtr alignment(openHalAlignment(inPath, &optionsParser, READ_ACCESS | WRITE_ACCESS));

    if (alignment->openGenome(root) == NULL) {
        cerr << "[halRemoveSubtree] Error: given root " << root << " not found in alignment" << endl;
        return 1;
    }
    
    if (!noMarkAncestors) {
        markAncestorsForUpdate(alignment, root);
    }
    
    // the main use case for this is prepping for halAppendSubtree, in which case we
    // want to leave the root node. 
    vector<string> children = alignment->getChildNames(root);
    if (children.empty()) {
        cerr << "[halRemoveSubtree] Warning: given root " << root << " is a leaf: doing nothing" << endl;
    }
    for (size_t i = 0; i < children.size(); ++i) {
        remove_recursive(alignment, children[i]);
    }
    return 0;
}
