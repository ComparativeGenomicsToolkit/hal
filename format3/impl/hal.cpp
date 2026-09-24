/*
 * hal: the HAL format 3 command line.  `hal <command> [options] <args>`.
 */
#include "hal3Cli.h"

#include <cstring>
#include <iostream>

using namespace std;

namespace {

struct Command {
    const char *name;
    const char *summary;
    int (*main)(int, char **);
};

const Command COMMANDS[] = {
#ifdef HAL3_CONVERT
    {"convert", "convert a HAL (HDF5 or mmap) alignment to format 3", hal3::main_convert},
#endif
    {"validate", "check an alignment's manifest, genome files and DNA", hal3::main_validate},
    {"liftover", "map BED intervals from one genome to another", hal3::main_liftover},
    {"maf", "export MAF blocks anchored on a reference genome", hal3::main_maf},
    {"fasta", "export genome sequences as FASTA", hal3::main_fasta},
};

void usage(ostream &os) {
    os << "usage: hal <command> [options] <args>\n\n"
          "HAL format 3: a directory holding a manifest, and per genome a .hgn genome file (its\n"
          "branch to the parent) and a .2bit of its DNA.\n\ncommands:\n";
    for (const Command &c : COMMANDS) {
        os << "  " << c.name;
        for (size_t k = strlen(c.name); k < 10; k++) os << ' ';
        os << c.summary << "\n";
    }
#ifndef HAL3_CONVERT
    os << "  (convert is not available: this hal was built without HDF5)\n";
#endif
    os << "\nhal <command> --help describes a command.\n";
}
}

int main(int argc, char **argv) {
    if (argc < 2 || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "help")) {
        usage(argc < 2 ? cerr : cout);
        return argc < 2 ? 1 : 0;
    }
    for (const Command &c : COMMANDS) {
        if (strcmp(argv[1], c.name)) continue;
        /* getopt reports errors under argv[0]: make that "hal <command>" */
        string prog = string("hal ") + c.name;
        argv[1] = &prog[0];
        return c.main(argc - 1, argv + 1);
    }
    if (!strcmp(argv[1], "convert"))
        cerr << "hal: convert is not available: this hal was built without HDF5" << endl;
    else
        cerr << "hal: unknown command '" << argv[1] << "'" << endl;
    usage(cerr);
    return 1;
}
