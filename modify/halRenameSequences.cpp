#include "hal.h"
#include "renameFile.h"

using namespace hal;
using namespace std;

static void initParser(CLParser &optionsParser) {
    optionsParser.setDescription("Rename the sequences of a genome in-place.");
    optionsParser.addArgument("halFile", "hal file");
    optionsParser.addArgument("genome", "genome to rename the sequences of");
    optionsParser.addArgument("renameFile", "Tab-separated file. First column: existing "
                                            "sequence name, second column: new sequence name."
                                            " Any sequences not provided will stay the same.");
}

int main(int argc, char *argv[]) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(optionsParser);
    string halPath, renamePath, genomeName;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        genomeName = optionsParser.getArgument<string>("genome");
        renamePath = optionsParser.getArgument<string>("renameFile");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        return 1;
    }

    AlignmentPtr alignment(openHalAlignment(halPath, &optionsParser, WRITE_ACCESS | READ_ACCESS));
    Genome *genome = alignment->openGenome(genomeName);
    if (genome == NULL) {
        throw hal_exception("Genome " + genomeName + " not found in alignment");
    }
    map<string, string> renameMap = ingestRenameFile(renamePath);

    for (map<string, string>::iterator it = renameMap.begin(); it != renameMap.end(); it++) {
        Sequence *sequence = genome->getSequenceCheck(it->first);

        sequence = genome->getSequence(it->second);
        if (sequence != NULL) {
            throw hal_exception("Attempting to rename sequence " + it->first + " to " + it->second + " failed: " + it->second +
                                " is already in the genome! Name it to something"
                                " temporary first");
        }
    }

    // Do the actual renaming now that we are relatively sure nothing
    // will go wrong.
    for (map<string, string>::iterator it = renameMap.begin(); it != renameMap.end(); it++) {
        cout << "Renaming " << it->first << " to " << it->second << endl;
        Sequence *sequence = genome->getSequence(it->first);
        sequence->setName(it->second);
    }

    alignment->close();
}
