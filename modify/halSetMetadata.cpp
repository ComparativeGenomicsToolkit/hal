// Simple utility to set metadata for a hal genome or alignment
#include "hal.h"

using namespace hal;
using namespace std;

static void initParser(CLParser &optionsParser) {
    optionsParser.setDescription("Set metadata for an alignment or genome");
    optionsParser.addArgument("halFile", "hal file to modify");
    optionsParser.addArgument("key", "metadata key");
    optionsParser.addArgument("value", "metadata value");
    optionsParser.addOption("genome", "genome to set metadata for instead of "
                                      "setting it for the entire alignment",
                            "");
}

int main(int argc, char *argv[]) {
    string halPath, key, value, genomeName;
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        key = optionsParser.getArgument<string>("key");
        value = optionsParser.getArgument<string>("value");
        genomeName = optionsParser.getOption<string>("genome");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        return 1;
    }

    AlignmentPtr alignment(openHalAlignment(halPath, &optionsParser,
                                            WRITE_ACCESS | READ_ACCESS));
    if (genomeName == "") {
        // No genome to set metadata for, so set the alignment-wide metadata.
        MetaData *metadata = alignment->getMetaData();
        metadata->set(key, value);
    } else {
        Genome *genome = alignment->openGenome(genomeName);
        MetaData *metadata = genome->getMetaData();
        metadata->set(key, value);
    }
    alignment->close();
    return 0;
}
