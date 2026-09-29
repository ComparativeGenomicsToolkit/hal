/*
 * hal merge: stitch small alignments (fragments -- e.g. one per Cactus subproblem) whose
 * trees meet in shared genomes into one alignment (see merge, hal3/hal3Writer.h).
 */
#include "hal3/hal3Writer.h"
#include "hal3Cli.h"
#include "hal3Manifest.h"
#include "hal3Store.h"

#include <getopt.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

namespace hal3 {

static const char *MERGE_USAGE =
    "usage: hal merge [options] <out.hal> <fragment>...\n"
    "Stitch alignments whose trees meet in shared genomes -- a fragment's root a leaf of\n"
    "another's, as a Cactus subproblem's ancestor is a leaf of its parent's -- into one.  A\n"
    "genome in several fragments takes its genome file from the one where it has a parent,\n"
    "its DNA from the one where it is the root; its DNA must agree between them.  Each file\n"
    "is copied once, checked against its fragment's manifest, and the result is validated;\n"
    "on any problem nothing is written.\n"
    "options:\n"
    "  --directory        write out.hal as a directory (files hard-linked from directory\n"
    "                     fragments where possible)\n"
    "  --meta KEY=VALUE   set alignment metadata (after the fragments'; repeatable)\n"
    "  --noDnaCheck       do not require a genome's DNA to agree between fragments\n"
    "  --checksums        also check every file's CRC-32 in the result\n"
    "  --tree NEWICK      the tree the result must be (same genomes and parents): its branch\n"
    "                     lengths and order of children, not the fragments'\n"
    "  -h, --help\n";

int main_merge(int argc, char **argv) {
    return runCommand("merge", [&]() {
        api::MergeOptions o;
        map<string, string> meta;
        enum { DIRECTORY = 1000, META, NO_DNA_CHECK, CHECKSUMS, TREE };
        static const option longOpts[] = {{"directory", no_argument, nullptr, DIRECTORY},
                                          {"meta", required_argument, nullptr, META},
                                          {"noDnaCheck", no_argument, nullptr, NO_DNA_CHECK},
                                          {"checksums", no_argument, nullptr, CHECKSUMS},
                                          {"tree", required_argument, nullptr, TREE},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case DIRECTORY: o.pack = false; break;
            case META: {
                string kv = optarg;
                size_t eq = kv.find('=');
                if (eq == string::npos || eq == 0) throw CliError("--meta needs KEY=VALUE, got " + kv);
                meta[kv.substr(0, eq)] = kv.substr(eq + 1);
                break;
            }
            case NO_DNA_CHECK: o.checkDna = false; break;
            case CHECKSUMS: o.checksums = true; break;
            case TREE: o.tree = optarg; break;
            case 'h': cout << MERGE_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind < 2) throw CliError("expected an output and at least one fragment");
        string out = argv[optind];
        vector<string> fragments(argv + optind + 1, argv + argc);
        api::merge(fragments, out, meta, o);
        Manifest m = Manifest::parse(Store::open(out)->read(MANIFEST_NAME), out);
        cerr << "hal merge: " << out << ": " << m.genomes.size() << " genomes from " << fragments.size() << " fragments" << endl;
        return 0;
    });
}
}
