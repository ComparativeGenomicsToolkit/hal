/*
 * hal update: change an alignment's genomes -- split, join, subset, replace, rename --
 * from alignments that exist (see hal3Update.h).
 */
#include "hal3Cli.h"
#include "hal3Update.h"

#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

namespace hal3 {

static const char *UPDATE_USAGE =
    "usage: hal update [options] <in.hal> <out.hal>\n"
    "       hal update --inPlace [options] <in.hal>\n"
    "Change an alignment's genomes, using alignments that exist: split, join, subset, replace,\n"
    "rename.  The operations apply in the order given; the result must be a valid alignment\n"
    "(one root; each genome's file mapping it to its parent's sequences) or nothing is\n"
    "written.  A genome's file is rewritten only when its name, its parent's name or its\n"
    "branch changes; the rest are copied untouched.  By default a new alignment is written\n"
    "and the input left alone.\n"
    "operations:\n"
    "  --remove G[,G...]      remove genomes.  A genome's children must go too, or be grafted\n"
    "                         or replaced below a parent that stays, in the same update\n"
    "  --prune G              remove G and everything below it\n"
    "  --extract G            keep only the subtree of G, which becomes the root (split)\n"
    "  --keep G[,G...]        keep only these genomes and the ancestors connecting them,\n"
    "                         rooted at their common ancestor (subset)\n"
    "  --graft DONOR          add alignment DONOR, whose root must be a genome here, below it\n"
    "                         (join)\n"
    "  --replace DONOR:G[,G...]  take these genomes' files (and parent and branch length)\n"
    "                         from alignment DONOR\n"
    "  --rename OLD=NEW       rename a genome\n"
    "  --branchLength G=X     set the length of the branch above G\n"
    "options:\n"
    "  --inPlace              change in.hal itself: for a single file, one crash-safe\n"
    "                         transaction whose cost is the new data, not the file\n"
    "  --directory            write out.hal as a directory\n"
    "  --dryRun               check the result and say what would change; write nothing\n"
    "  -h, --help\n";

namespace {
vector<string> commas(const string &s) {
    vector<string> v;
    stringstream ss(s);
    string t;
    while (getline(ss, t, ','))
        if (!t.empty()) v.push_back(t);
    return v;
}
pair<string, string> split(const string &s, char sep, const char *what) {
    size_t k = s.rfind(sep);
    if (k == string::npos || k == 0 || k + 1 == s.size()) throw CliError(string(what) + " needs A" + sep + "B, not " + s);
    return {s.substr(0, k), s.substr(k + 1)};
}
}

int main_update(int argc, char **argv) {
    return runCommand("update", [&]() {
        enum { REMOVE = 1000, PRUNE, EXTRACT, KEEP, GRAFT, REPLACE, RENAME, BRANCH_LENGTH, IN_PLACE, DIRECTORY, DRY_RUN };
        static const option longOpts[] = {{"remove", required_argument, nullptr, REMOVE},
                                          {"prune", required_argument, nullptr, PRUNE},
                                          {"extract", required_argument, nullptr, EXTRACT},
                                          {"keep", required_argument, nullptr, KEEP},
                                          {"graft", required_argument, nullptr, GRAFT},
                                          {"replace", required_argument, nullptr, REPLACE},
                                          {"rename", required_argument, nullptr, RENAME},
                                          {"branchLength", required_argument, nullptr, BRANCH_LENGTH},
                                          {"inPlace", no_argument, nullptr, IN_PLACE},
                                          {"directory", no_argument, nullptr, DIRECTORY},
                                          {"dryRun", no_argument, nullptr, DRY_RUN},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        vector<pair<int, string>> ops;
        bool inPlace = false, directory = false, dryRun = false;
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case IN_PLACE: inPlace = true; break;
            case DIRECTORY: directory = true; break;
            case DRY_RUN: dryRun = true; break;
            case 'h': cout << UPDATE_USAGE; return 0;
            default:
                if (c < REMOVE || c > BRANCH_LENGTH) throw CliError("bad option");
                ops.push_back({c, optarg});
            }
        }
        int want = inPlace || dryRun ? 1 : 2;
        if (argc - optind != want && !(dryRun && argc - optind == 2))
            throw CliError("expected " + to_string(want) + " argument" + (want > 1 ? "s" : "") + ", got " + to_string(argc - optind));
        if (inPlace && directory) throw CliError("--inPlace changes the input in its own form; --directory is for new output");
        if (ops.empty()) throw CliError("no operations given");
        Update u(argv[optind]);
        for (auto &op : ops) {
            switch (op.first) {
            case REMOVE: u.remove(commas(op.second)); break;
            case PRUNE: u.prune(op.second); break;
            case EXTRACT: u.extract(op.second); break;
            case KEEP: u.keep(commas(op.second)); break;
            case GRAFT: u.graft(op.second); break;
            case REPLACE: {
                auto dg = split(op.second, ':', "--replace");
                u.replace(dg.first, commas(dg.second));
                break;
            }
            case RENAME: {
                auto ft = split(op.second, '=', "--rename");
                u.rename(ft.first, ft.second);
                break;
            }
            case BRANCH_LENGTH: {
                auto gl = split(op.second, '=', "--branchLength");
                u.setBranchLength(gl.first, parseDouble("branchLength", gl.second.c_str()));
                break;
            }
            }
        }
        if (dryRun) {
            vector<string> p = u.check();
            for (const string &l : u.log()) cerr << "hal update: " << l << endl;
            if (p.empty()) {
                cerr << "hal update: the result is a valid alignment (nothing written)" << endl;
                return 0;
            }
            for (const string &x : p) cerr << "hal update: problem: " << x << endl;
            return 1;
        }
        if (inPlace) u.writeInPlace();
        else u.write(argv[optind + 1], directory);
        for (const string &l : u.log()) cerr << "hal update: " << l << endl;
        return 0;
    });
}
}
