/*
 * hal3Liftover: halLiftover for HAL format 3.  Same arguments and the same BED output;
 * no HAL API and no HDF5.
 */
#include "hal3Alignment.h"
#include "hal3Liftover.h"

#include <cstring>
#include <fstream>
#include <iostream>

using namespace std;

static void usage() {
    cerr << "usage: hal3Liftover alignment_dir srcGenome srcBed tgtGenome tgtBed [--noDupes] [--append] [--bedType N]\n"
            "Map BED intervals from srcGenome to tgtGenome, as halLiftover does.  srcBed may be\n"
            "stdin and tgtBed stdout.\n"
            "  --noDupes     do not map between duplications: one target copy per base\n"
            "  --append      append to tgtBed\n"
            "  --bedType N   number of standard BED columns (3 to 12); more are passed through.\n"
            "                By default taken from the first line.\n";
    exit(1);
}

int main(int argc, char **argv) {
    vector<string> args;
    bool noDupes = false, append = false;
    int bedType = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--noDupes")) noDupes = true;
        else if (!strcmp(argv[i], "--append")) append = true;
        else if (!strcmp(argv[i], "--bedType") && i + 1 < argc) {
            bedType = atoi(argv[++i]);
            if (bedType < 3 || bedType > 12) {
                cerr << "--bedType must be between 3 and 12" << endl;
                return 1;
            }
        } else if (!strcmp(argv[i], "--outPSL") || !strcmp(argv[i], "--outPSLWithName") || !strcmp(argv[i], "--coalescenceLimit")) {
            cerr << "hal3Liftover: " << argv[i] << " is not implemented yet" << endl;
            return 1;
        } else if (argv[i][0] == '-' && argv[i][1] == '-') usage();
        else args.push_back(argv[i]);
    }
    if (args.size() != 5) usage();
    try {
        hal3::Alignment alignment(args[0]);
        const hal3::Genome *src = alignment.genome(args[1]), *tgt = alignment.genome(args[3]);
        if (!src) throw runtime_error("srcGenome, " + args[1] + ", not found in alignment");
        if (!tgt) throw runtime_error("tgtGenome, " + args[3] + ", not found in alignment");
        ifstream inFile;
        istream *in = &cin;
        if (args[2] != "stdin") {
            inFile.open(args[2]);
            if (!inFile) throw runtime_error("Error opening srcBed, " + args[2]);
            in = &inFile;
        }
        ofstream outFile;
        ostream *out = &cout;
        if (args[4] != "stdout") {
            outFile.open(args[4], append ? ios::out | ios::app : ios::out);
            if (!outFile) throw runtime_error("Error opening tgtBed, " + args[4]);
            out = &outFile;
        }
        hal3::BedLiftover liftover(src, tgt, !noDupes, out);
        liftover.convert(in, bedType);
    } catch (exception &e) {
        cerr << "hal3Liftover: " << e.what() << endl;
        return 1;
    }
    return 0;
}
