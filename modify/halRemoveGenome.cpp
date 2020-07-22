#include "hal.h"
#include "halCLParser.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("inFile", "existing tree");
    optionsParser.addArgument("deleteNode", "(leaf) genome to delete");
    optionsParser.addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                                                   " update",
                                false);
}

int main(int argc, char *argv[]) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    string inPath, deleteNode;
    bool noMarkAncestors;
    try {
        optionsParser.parseOptions(argc, argv);
        inPath = optionsParser.getArgument<string>("inFile");
        deleteNode = optionsParser.getArgument<string>("deleteNode");
        noMarkAncestors = optionsParser.getFlag("noMarkAncestors");
    } catch (exception &e) {
        optionsParser.printUsage(cerr);
        return 1;
    }
    AlignmentPtr alignment(openHalAlignment(inPath, &optionsParser, READ_ACCESS | WRITE_ACCESS));
    if (!noMarkAncestors) {
        markAncestorsForUpdate(alignment, deleteNode);
    }
    alignment->removeGenome(deleteNode);
    return 0;
}
