#include "hal.h"
#include "markAncestors.h"

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("inFile", "existing tree");
    optionsParser.addArgument("botAlignmentFile", "tree containing insert, its "
                                                  "proper bottom segments, and the new leaf genome");
    optionsParser.addArgument("topAlignmentFile", "tree containing insert, its "
                                                  "parent, and its proper top segments");
    optionsParser.addArgument("parentName", "insert's future parent");
    optionsParser.addArgument("insertName", "insert name");
    optionsParser.addArgument("childName", "insert's future child");
    optionsParser.addArgument("leafName", "name of new leaf genome");
    optionsParser.addArgument("upperBranchLength", "length of branch from parent"
                                                   " to insert");
    optionsParser.addArgument("leafBranchLength", "leaf branch length");
    optionsParser.addOptionFlag("noMarkAncestors", "don't mark ancestors for"
                                                   " update",
                                false);
}

int main(int argc, char *argv[]) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    string inPath, botAlignmentPath, topAlignmentPath, parentName, insertName, childName, leafName;
    double upperBranchLength, leafBranchLength;
    bool noMarkAncestors;
    try {
        optionsParser.parseOptions(argc, argv);
        inPath = optionsParser.getArgument<string>("inFile");
        botAlignmentPath = optionsParser.getArgument<string>("botAlignmentFile");
        topAlignmentPath = optionsParser.getArgument<string>("topAlignmentFile");
        parentName = optionsParser.getArgument<string>("parentName");
        insertName = optionsParser.getArgument<string>("insertName");
        childName = optionsParser.getArgument<string>("childName");
        leafName = optionsParser.getArgument<string>("leafName");
        upperBranchLength = optionsParser.getArgument<double>("upperBranchLength");
        leafBranchLength = optionsParser.getArgument<double>("leafBranchLength");
        noMarkAncestors = optionsParser.getFlag("noMarkAncestors");
    } catch (exception &e) {
        optionsParser.printUsage(cerr);
        return 1;
    }
    AlignmentPtr mainAlignment(openHalAlignment(inPath, &optionsParser, READ_ACCESS | WRITE_ACCESS));
    AlignmentConstPtr botAlignment(openHalAlignment(botAlignmentPath, &optionsParser));
    AlignmentConstPtr topAlignment(openHalAlignment(topAlignmentPath, &optionsParser));

    // Add in the insert genome.
    mainAlignment->insertGenome(insertName, parentName, childName, upperBranchLength);

    // Add in the new leaf genome and its dimensions (so that it has the proper number of sequences).
    mainAlignment->addLeafGenome(leafName, insertName, leafBranchLength);
    Genome *outLeafGenome = mainAlignment->openGenomeCheck(leafName);
    const Genome *inLeafGenome = botAlignment->openGenomeCheck(leafName);
    inLeafGenome->copyDimensions(outLeafGenome);

    // Insert the new intermediate node.
    Genome *insertGenome = mainAlignment->openGenomeCheck(insertName);
    const Genome *topInsertGenome = topAlignment->openGenomeCheck(insertName);
    const Genome *botInsertGenome = botAlignment->openGenomeCheck(insertName);
    topInsertGenome->copyDimensions(insertGenome);
    topInsertGenome->copyTopDimensions(insertGenome);
    botInsertGenome->copyBottomDimensions(insertGenome);
    topInsertGenome->copySequence(insertGenome);
    topInsertGenome->copyTopSegments(insertGenome);
    topInsertGenome->copyMetadata(insertGenome);
    botInsertGenome->copyBottomSegments(insertGenome);
    insertGenome->fixParseInfo();

    // Copy the bottom segments for the parent genome from the top alignment.
    Genome *parentGenome = mainAlignment->openGenomeCheck(parentName);
    const Genome *botParentGenome = topAlignment->openGenomeCheck(parentName);
    botParentGenome->copyBottomDimensions(parentGenome);
    botParentGenome->copyBottomSegments(parentGenome);
    parentGenome->fixParseInfo();

    // Fix the parent's other children as well.
    vector<string> allChildren = mainAlignment->getChildNames(parentName);
    for (size_t i = 0; i < allChildren.size(); i++) {
        if (allChildren[i] != insertName) {
            Genome *outGenome = mainAlignment->openGenomeCheck(allChildren[i]);
            const Genome *topSegmentsGenome = topAlignment->openGenomeCheck(allChildren[i]);
            topSegmentsGenome->copyTopDimensions(outGenome);
            topSegmentsGenome->copyTopSegments(outGenome);
            outGenome->fixParseInfo();
        }
    }

    // Copy the top segments for the child genome from the bottom alignment.
    Genome *childGenome = mainAlignment->openGenomeCheck(childName);
    const Genome *topChildGenome = botAlignment->openGenomeCheck(childName);
    topChildGenome->copyTopDimensions(childGenome);
    topChildGenome->copyTopSegments(childGenome);
    childGenome->fixParseInfo();

    // Copy the entire genome for the leaf from the bottom alignment.
    inLeafGenome->copy(outLeafGenome);
    if (!noMarkAncestors) {
        markAncestorsForUpdate(mainAlignment, insertName);
    }
}
