#include "hal.h"
#include "halAlignmentInstance.h"
#include "halCLParser.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("mainFile", "destination tree");
    optionsParser.addArgument("appendFile", "alignment containing the tree to be"
                                            " appended");
    optionsParser.addArgument("parentName", "node to be added to");
    optionsParser.addArgument("rootName", "name of subtree root");
    optionsParser.addOption("bridgeFile", "alignment containing parent,"
                                          " subtree root, and its future siblings, if any "
                                          "(required if not merging appended and appendee "
                                          "nodes)",
                            "");
    optionsParser.addOption("branchLength", "branch length between appended and "
                                            "appendee nodes",
                            0.0);
    optionsParser.addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                                                   " update",
                                false);
    optionsParser.addOptionFlag("merge", "merge appended root and node that is appended to", false);
}

static void closeWithNeighbours(AlignmentConstPtr alignment, const Genome* genome) {
    // Genome::copy() can open parent and child genomes
    // but we want to make absolutely sure they're closed afterwards here so we can
    // run with --inMemory
    const Genome* parent = genome->getParent();
    if (parent) {
        alignment->closeGenome(parent);
    }
    size_t ccount = genome->getNumChildren();
    for (size_t i = 0; i < ccount; ++i) {
        const Genome* child = genome->getChild(i);
        alignment->closeGenome(child);
    }
    alignment->closeGenome(genome);
}

void addSubtree(AlignmentPtr mainAlignment, AlignmentConstPtr appendAlignment, string currNode) {
    vector<string> children = appendAlignment->getChildNames(currNode);
    for (size_t i = 0; i < children.size(); i++) {
        Genome* testGenome = mainAlignment->openGenome(children[i]);
        if (testGenome != NULL) {
            // Special case for a node that is being merged; we don't want
            // to append children that already exist.
            mainAlignment->closeGenome(testGenome);
            continue;
        }
        Genome *mainChildGenome =
            mainAlignment->addLeafGenome(children[i], currNode, appendAlignment->getBranchLength(currNode, children[i]));
        const Genome *appendChildGenome = appendAlignment->openGenome(children[i]);
        cerr << "[halAppendSubtree] Copying " << children[i] << endl;
        appendChildGenome->copy(mainChildGenome);
        closeWithNeighbours(mainAlignment, mainChildGenome);        
        closeWithNeighbours(appendAlignment, appendChildGenome);        
        addSubtree(mainAlignment, appendAlignment, children[i]);
    }
    Genome *outGenome = mainAlignment->openGenome(currNode);
    const Genome *inGenome = appendAlignment->openGenome(currNode);
    inGenome->copyBottomDimensions(outGenome);
    inGenome->copyBottomSegments(outGenome);
    outGenome->fixParseInfo();
    closeWithNeighbours(mainAlignment, outGenome);
    closeWithNeighbours(appendAlignment, inGenome);    
}

int main(int argc, char *argv[]) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    string mainPath, appendPath, bridgePath, parentName, rootName;
    double branchLength;
    bool noMarkAncestors;
    bool merge;
    bool inMemory;
    try {
        optionsParser.parseOptions(argc, argv);
        mainPath = optionsParser.getArgument<string>("mainFile");
        appendPath = optionsParser.getArgument<string>("appendFile");
        bridgePath = optionsParser.getOption<string>("bridgeFile");
        parentName = optionsParser.getArgument<string>("parentName");
        rootName = optionsParser.getArgument<string>("rootName");
        branchLength = optionsParser.getOption<double>("branchLength");
        noMarkAncestors = optionsParser.getFlag("noMarkAncestors");
        merge = optionsParser.getFlag("merge");
        inMemory = optionsParser.getFlag("hdf5InMemory") || optionsParser.getFlag("inMemory");
    } catch (exception &e) {
        optionsParser.printUsage(cerr);
        return 1;
    }
    AlignmentPtr mainAlignment(openHalAlignment(mainPath, &optionsParser, READ_ACCESS | WRITE_ACCESS));
    AlignmentConstPtr appendAlignment(openHalAlignment(appendPath, &optionsParser));
    AlignmentConstPtr bridgeAlignment;

    if (!inMemory) {
        cerr << "[halAppendSubtree] Warning this tool requires --hdf5InMemory to be set in order to run in a reasonable time" << endl;
    }

    if (!merge) {
        if (bridgePath == "") {
            throw hal_exception("need a bridge alignment if not merging nodes");
        }
        bridgeAlignment = AlignmentConstPtr(openHalAlignment(bridgePath, &optionsParser));
        Genome *mainAppendedRoot = mainAlignment->addLeafGenome(rootName, parentName, branchLength);
        const Genome *appendAppendedRoot = appendAlignment->openGenome(rootName);
        cerr << "[halAppendSubtree] Copying " << rootName << endl;
        appendAppendedRoot->copy(mainAppendedRoot);
        appendAlignment->closeGenome(appendAppendedRoot);
    } else {
        // the bridge alignment is equivalent to the append alignment in this case
        // (the append alignment will contain at least all the information that
        // the bridge alignment would)
        bridgeAlignment = appendAlignment;
        if (parentName != rootName) {
            throw hal_exception("parent name must be equal to root name if "
                                "--merge option is given");
        }
        assert(branchLength == 0.0);
    }
    addSubtree(mainAlignment, appendAlignment, rootName);
    
    // Need proper bottom segments for parent genome
    // 1) Copy dimensions for the genome that's being appended
    Genome *mainParentGenome = mainAlignment->openGenome(parentName);
    const Genome *bridgeParentGenome = bridgeAlignment->openGenome(parentName);
    bridgeParentGenome->copyBottomDimensions(mainParentGenome);
    
    // 2) Copy top dimensions for the children
    vector<string> children = bridgeAlignment->getChildNames(parentName);
    for (size_t i = 0; i < children.size(); i++) {
        Genome *mainChildGenome = mainAlignment->openGenome(children[i]);
        const Genome *bridgeChildGenome = bridgeAlignment->openGenome(children[i]);
        bridgeChildGenome->copyTopDimensions(mainChildGenome);
    }

    // 3) copy bottom segments for the genome that's being appended
    bridgeParentGenome->copyBottomSegments(mainParentGenome);
    mainParentGenome->fixParseInfo();
    
    // 4) copu top segments for its children
    for (size_t i = 0; i < children.size(); i++) {
        Genome *mainChildGenome = mainAlignment->openGenome(children[i]);
        const Genome *bridgeChildGenome = bridgeAlignment->openGenome(children[i]);
        bridgeChildGenome->copyTopSegments(mainChildGenome);
        mainChildGenome->fixParseInfo();
        mainAlignment->closeGenome(mainChildGenome);
        bridgeAlignment->closeGenome(bridgeChildGenome);
    } 
    
    mainAlignment->closeGenome(mainParentGenome);
    bridgeAlignment->closeGenome(bridgeParentGenome);


    if (!noMarkAncestors) {
        markAncestorsForUpdate(mainAlignment, rootName);
    }
    return 0;
}
