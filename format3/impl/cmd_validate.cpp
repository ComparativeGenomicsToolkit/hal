/*
 * hal validate: check a HAL format 3 alignment, pack or directory (the checks are
 * api::validate's, hal3/hal3Writer.h).
 *
 * Uses only the manifest code, liftmap and twobit64 -- no HAL API, no HDF5.
 *
 *   manifest    parses; the Newick tree and the roster name the same genomes with the
 *               same parents and branch lengths; the root is first
 *   files       every .hgn and .2bit exists with the size and CRC-32 the manifest
 *               records (skip the CRC with --noCrc)
 *   .hgn        profile hal.branch; hal.genome and hal.parent match the manifest;
 *               liftmap's full integrity check (every checksum, every chunk decoded and
 *               inside its directory ranges); axis a is a function, i.e. every base has
 *               at most one parent position; the root has no runs and no axis b
 *   tables      axis a matches the manifest's length and sequence count; axis b is the
 *               parent's axis a, name for name and length for length; the .2bit holds
 *               the same sequences as axis a
 *   members     each exists with its size and CRC-32; a likelihood (lik) member is for
 *               the genome's own sequences, and every chunk of it is sound and decodes
 *   orphans     no member the manifest does not name
 *   pack        every member's bytes match the CRC-32 in the pack's index
 */
#include "hal3/hal3Writer.h"
#include "hal3Cli.h"
#include "hal3Manifest.h"
#include "hal3Store.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace hal3 {

static const char *VALIDATE_USAGE =
    "usage: hal validate [options] <alignment>\n"
    "Check an alignment: manifest, tree, files, checksums, and every genome file against its\n"
    "parent's and its DNA.\n"
    "options:\n"
    "  --noCrc      skip the whole-file checksums\n"
    "  -h, --help\n";

static int validate(int argc, char **argv) {
    bool doCrc = true;
    enum { NO_CRC = 1000 };
    static const option longOpts[] = {{"noCrc", no_argument, nullptr, NO_CRC},
                                      {"help", no_argument, nullptr, 'h'},
                                      {nullptr, 0, nullptr, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
        switch (c) {
        case NO_CRC: doCrc = false; break;
        case 'h': cout << VALIDATE_USAGE; return 0;
        default: throw CliError("bad option");
        }
    }
    if (argc - optind != 1) throw CliError("expected 1 argument, got " + to_string(argc - optind));
    string dir = argv[optind];
    Manifest m;
    try {
        m = Manifest::parse(Store::open(dir)->read(MANIFEST_NAME), dir);
    } catch (exception &e) {
        cerr << "hal validate: " << e.what() << endl;
        return 1;
    }
    api::ValidateOptions vo;
    vo.checksums = doCrc;
    vector<string> problems = api::validate(dir, vo);
    for (size_t i = 0; i < problems.size() && i < 50; i++) cerr << "  FAIL " << problems[i] << endl;
    cout << dir << ": " << m.genomes.size() << " genomes, " << (problems.empty() ? "ok" : "FAILED (" + to_string(problems.size()) + " problems)")
         << endl;
    return problems.empty() ? 0 : 1;
}

int main_validate(int argc, char **argv) { return runCommand("validate", [&]() { return validate(argc, argv); }); }
}
