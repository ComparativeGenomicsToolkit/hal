/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halStats.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace hal;

int main(int argc, char **argv) {
    CLParser optionsParser;
    optionsParser.addArgument("halFile", "path to hal file to validate");
    optionsParser.addOption("genome", "specific genome to validate instead of entire file", "");
    optionsParser.setDescription("Check if hal database is valid");
    string path, genomeName;
    try {
        optionsParser.parseOptions(argc, argv);
        path = optionsParser.getArgument<string>("halFile");
        genomeName = optionsParser.getOption<string>("genome");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
    try {
        AlignmentConstPtr alignment(openHalAlignment(path, &optionsParser));
        if (genomeName == "") {
            validateAlignment(alignment.get());
        } else {
            const Genome *genome = alignment->openGenome(genomeName);
            validateGenome(genome);
        }
    } catch (hal_exception &e) {
        cerr << "hal exception caught: " << e.what() << endl;
        return 1;
    } catch (exception &e) {
        cerr << "Exception caught: " << e.what() << endl;
        return 1;
    }
    cout << "\nFile valid" << endl;

    return 0;
}
