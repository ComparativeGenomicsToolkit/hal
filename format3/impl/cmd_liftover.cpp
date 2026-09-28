/*
 * hal liftover: map BED intervals between genomes, as halLiftover (same positional
 * arguments, same BED output).
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3BedLiftover.h"

#include <fstream>
#include <getopt.h>
#include <iostream>

using namespace std;

namespace hal3 {

static const char *LIFTOVER_USAGE =
    "usage: hal liftover [options] <alignment> <srcGenome> <srcBed> <tgtGenome> <tgtBed>\n"
    "Map BED intervals from srcGenome to tgtGenome, as halLiftover does.  srcBed may be stdin,\n"
    "tgtBed stdout.  Output is halLiftover's; its lines may be split differently where\n"
    "paralogous copies overlap on the target, which depended on HAL's segmentation.\n"
    "options:\n"
    "  --noDupes        do not map between duplications: one target copy per base\n"
    "  --append         append to tgtBed\n"
    "  --bedType N      number of standard BED columns (3 to 12); more are passed through.\n"
    "                   By default taken from the first line.\n"
    "  -h, --help\n";

int main_liftover(int argc, char **argv) {
    return runCommand("liftover", [&]() {
        bool noDupes = false, append = false;
        int bedType = 0;
        enum { NO_DUPES = 1000, APPEND, BED_TYPE };
        static const option longOpts[] = {{"noDupes", no_argument, nullptr, NO_DUPES},
                                          {"append", no_argument, nullptr, APPEND},
                                          {"bedType", required_argument, nullptr, BED_TYPE},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case NO_DUPES: noDupes = true; break;
            case APPEND: append = true; break;
            case BED_TYPE:
                bedType = (int)parseInt("bedType", optarg);
                if (bedType < 3 || bedType > 12) throw CliError("--bedType must be between 3 and 12");
                break;
            case 'h': cout << LIFTOVER_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 5) throw CliError("expected 5 arguments, got " + to_string(argc - optind));
        string dir = argv[optind], srcName = argv[optind + 1], srcBed = argv[optind + 2], tgtName = argv[optind + 3],
               tgtBed = argv[optind + 4];
        Alignment alignment(dir);
        const Genome *src = alignment.genome(srcName), *tgt = alignment.genome(tgtName);
        if (!src) throw runtime_error("srcGenome, " + srcName + ", not found in alignment");
        if (!tgt) throw runtime_error("tgtGenome, " + tgtName + ", not found in alignment");
        ifstream inFile;
        istream *in = &cin;
        if (srcBed != "stdin") {
            inFile.open(srcBed);
            if (!inFile) throw runtime_error("Error opening srcBed, " + srcBed);
            in = &inFile;
        }
        ofstream outFile;
        ostream *out = &cout;
        if (tgtBed != "stdout") {
            outFile.open(tgtBed, append ? ios::out | ios::app : ios::out);
            if (!outFile) throw runtime_error("Error opening tgtBed, " + tgtBed);
            out = &outFile;
        }
        BedLiftover liftover(src, tgt, !noDupes, out);
        liftover.convert(in, bedType);
        out->flush();
        if (!*out) throw runtime_error("write to " + tgtBed + " failed");
        return 0;
    });
}
}
