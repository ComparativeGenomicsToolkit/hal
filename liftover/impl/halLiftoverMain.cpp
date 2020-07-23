/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halBlockLiftover.h"
#include "halColumnLiftover.h"
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;
using namespace hal;

static void initParser(CLParser &optionsParser) {
    optionsParser.addArgument("halFile", "input hal file");
    optionsParser.addArgument("srcGenome", "source genome name");
    optionsParser.addArgument("srcBed", "path of input bed file.  set as stdin "
                                        "to stream from standard input");
    optionsParser.addArgument("tgtGenome", "target genome name");
    optionsParser.addArgument("tgtBed", "path of output bed file.  set as stdout"
                                        " to stream to standard output.");
    optionsParser.addOptionFlag("noDupes", "do not map between duplications in"
                                           " graph.",
                                false);
    optionsParser.addOptionFlag("append", "append results to tgtBed", false);
    optionsParser.addOption("coalescenceLimit", "coalescence limit genome:"
                                                " the genome at or above the MRCA of source"
                                                " and target at which we stop looking for"
                                                " homologies (default: MRCA)",
                            "");
    optionsParser.addOptionFlag("outPSL", "write output in PSL instead of bed format",
                                false);
    optionsParser.addOptionFlag("outPSLWithName", "write output as input BED name followed by PSL line instead of "
                                                  "bed format",
                                false);
    optionsParser.addOption("bedType", "number of standard columns (3 to 12), columns beyond this are passed "
                            "through.  This only needs to be specified for BEDs with less than 12 columns and "
                            "having non-standard extra columns.", 0);
    optionsParser.setDescription("Map BED or PSL genome interval coordinates between "
                                 "two genomes.");
}

int main(int argc, char **argv) {
    CLParser optionsParser;
    initParser(optionsParser);

    string halPath;
    string srcGenomeName;
    string srcBedPath;
    string tgtGenomeName;
    string tgtBedPath;
    string coalescenceLimitName;
    bool noDupes;
    bool append;
    int bedType;
    bool outPSL;
    bool outPSLWithName;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        srcGenomeName = optionsParser.getArgument<string>("srcGenome");
        srcBedPath = optionsParser.getArgument<string>("srcBed");
        tgtGenomeName = optionsParser.getArgument<string>("tgtGenome");
        tgtBedPath = optionsParser.getArgument<string>("tgtBed");
        coalescenceLimitName = optionsParser.getOption<string>("coalescenceLimit");
        noDupes = optionsParser.getFlag("noDupes");
        append = optionsParser.getFlag("append");
        if (optionsParser.specifiedOption("bedType")) {
            bedType = optionsParser.getOption<int>("bedType");
            if ((bedType < 3) or (bedType > 12)) {
                throw hal_exception("--bedType must be between 3 and 12");
            }
        } else {
            bedType = 0;
        }
        outPSL = optionsParser.getFlag("outPSL");
        outPSLWithName = optionsParser.getFlag("outPSLWithName");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    try {
        if (outPSLWithName == true) {
            outPSL = true;
        }
        AlignmentConstPtr alignment(openHalAlignment(halPath, &optionsParser));
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("hal alignment is empty");
        }

        const Genome *srcGenome = alignment->openGenome(srcGenomeName);
        const Genome *tgtGenome = alignment->openGenome(tgtGenomeName);

        const Genome *coalescenceLimit = NULL;
        if (coalescenceLimitName != "") {
            coalescenceLimit = alignment->openGenome(coalescenceLimitName);
        }

        ifstream srcBed;
        istream *srcBedPtr;
        if (srcBedPath == "stdin") {
            srcBedPtr = &cin;
        } else {
            srcBed.open(srcBedPath.c_str());
            srcBedPtr = &srcBed;
            if (!srcBed) {
                throw hal_exception("Error opening srcBed, " + srcBedPath);
            }
        }

        ios_base::openmode mode = append ? ios::out | ios::app : ios_base::out;
        ofstream tgtBed;
        ostream *tgtBedPtr;
        if (tgtBedPath == "stdout") {
            tgtBedPtr = &cout;
        } else {
            tgtBed.open(tgtBedPath.c_str(), mode);
            tgtBedPtr = &tgtBed;
            if (!tgtBed) {
                throw hal_exception("Error opening tgtBed, " + tgtBedPath);
            }
        }

        BlockLiftover liftover;
        liftover.convert(alignment, srcGenome, srcBedPtr, tgtGenome, tgtBedPtr, bedType,
                         !noDupes, outPSL, outPSLWithName, coalescenceLimit);


    } catch (hal_exception &e) {
        cerr << "hal exception caught: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }

    return 0;
}
