/*
 * hal import-c2h: a Cactus subproblem in the .c2h form to a fragment (see
 * hal3/hal3Cactus.h), for hal merge to stitch with its neighbours'.
 */
#include "hal3/hal3Cactus.h"
#include "hal3Cli.h"
#include "hal3Manifest.h"
#include "hal3Store.h"

#include <getopt.h>
#include <iostream>
#include <string>

using namespace std;

namespace hal3 {

static const char *IMPORT_USAGE =
    "usage: hal import-c2h [options] <in.c2h> <in.fa> <tree> <fragment>\n"
    "Write a Cactus subproblem -- cactus_consolidated's .c2h, its FASTA, and the subproblem's\n"
    "tree (the ancestor and its children, a Newick string or file, as for\n"
    "halAppendCactusSubtree) -- as a fragment directory: the ancestor its root, each child with\n"
    "its branch.  Genomes the tree does not name (outgroups) are left out.  hal merge stitches\n"
    "the fragments of a whole alignment together.\n"
    "options:\n"
    "  --presorted   runs arrive in genome order (as Cactus writes them): stream, do not sort\n"
    "  -h, --help\n";

int main_import(int argc, char **argv) {
    return runCommand("import-c2h", [&]() {
        api::GenomeWriterOptions o;
        o.source = "cactus";
        enum { PRESORTED = 1000 };
        static const option longOpts[] = {{"presorted", no_argument, nullptr, PRESORTED},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case PRESORTED: o.presorted = true; break;
            case 'h': cout << IMPORT_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 4) throw CliError("expected 4 arguments, got " + to_string(argc - optind));
        string out = argv[optind + 3];
        api::importC2h(argv[optind], argv[optind + 1], argv[optind + 2], out, o);
        Manifest m = Manifest::parse(Store::open(out)->read(MANIFEST_NAME), out);
        cerr << "hal import-c2h: " << out << ": " << m.genomes.size() << " genomes, root " << m.genomes[0].name << endl;
        return 0;
    });
}
}
