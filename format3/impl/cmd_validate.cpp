/*
 * hal validate: check a HAL format 3 alignment, pack or directory.
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
#include "hal3Cli.h"
#include "hal3Likelihood.h"
#include "hal3Manifest.h"
#include "hal3Store.h"
#include "liftmap.h"
#include "twobit64.h"

#include <cmath>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
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

static int errors = 0;
static void fail(const string &what) {
    if (errors++ < 50) cerr << "  FAIL " << what << endl;
}

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
    errors = 0;

    hal3::Manifest m;
    unique_ptr<hal3::Store> store;
    try {
        store = hal3::Store::open(dir);
        m = hal3::Manifest::parse(store->read(hal3::MANIFEST_NAME), dir);
    } catch (exception &e) {
        cerr << "hal validate: " << e.what() << endl;
        return 1;
    }

    /* ---- the tree against the roster */
    try {
        vector<hal3::NewickNode> nodes = hal3::parseNewick(m.newick);
        set<string> seen;
        for (const hal3::NewickNode &n : nodes) {
            if (n.name.empty()) { fail("the tree has an unnamed node"); continue; }
            if (!seen.insert(n.name).second) fail("the tree names " + n.name + " twice");
            const hal3::GenomeEntry *g = m.find(n.name);
            if (!g) { fail("genome " + n.name + " is in the tree but not the roster"); continue; }
            string parent = n.parent < 0 ? "" : nodes[n.parent].name;
            if (g->parent != parent) fail(n.name + ": parent is " + parent + " in the tree, " + g->parent + " in the roster");
            if (n.parent >= 0 && n.branchLength >= 0 && fabs(n.branchLength - g->branchLength) > 1e-9 * max(1.0, fabs(n.branchLength)))
                fail(n.name + ": branch length differs between tree and roster");
        }
        for (const hal3::GenomeEntry &g : m.genomes)
            if (!seen.count(g.name)) fail("genome " + g.name + " is in the roster but not the tree");
        if (!m.genomes.empty() && !m.genomes[0].parent.empty()) fail("the first genome listed is not the root");
    } catch (exception &e) {
        fail(e.what());
    }

    /* ---- each genome */
    map<string, vector<pair<string, uint64_t>>> tables;   /* genome -> its axis a */
    set<string> named;
    for (const hal3::GenomeEntry &g : m.genomes) {
        string gname = "genome " + g.name;
        named.insert(g.hgnPath);
        if (!g.dnaPath.empty()) named.insert(g.dnaPath);
        for (const hal3::MemberRef &mr : g.members) {
            if (!named.insert(mr.path).second) fail(gname + ": its " + mr.role + " member " + mr.path + " is named twice");
            if (!store->has(mr.path)) {
                fail(gname + ": its " + mr.role + " member " + mr.path + " is missing");
                continue;
            }
            if (doCrc) {
                try {
                    if ((uint64_t)store->size(mr.path) != mr.bytes || store->crc(mr.path) != mr.crc)
                        fail(gname + ": " + mr.path + " size or CRC differs from the manifest");
                } catch (exception &e) {
                    fail(gname + ": " + e.what());
                }
            }
        }
        if (doCrc) {
            try {
                if ((uint64_t)store->size(g.hgnPath) != g.hgnBytes || store->crc(g.hgnPath) != g.hgnCrc)
                    fail(gname + ": " + g.hgnPath + " size or CRC differs from the manifest");
                if (!g.dnaPath.empty() && ((uint64_t)store->size(g.dnaPath) != g.dnaBytes || store->crc(g.dnaPath) != g.dnaCrc))
                    fail(gname + ": " + g.dnaPath + " size or CRC differs from the manifest");
            } catch (exception &e) {
                fail(gname + ": " + e.what());
            }
        }
        if (!store->has(g.hgnPath)) { fail(gname + ": " + g.hgnPath + " is missing"); continue; }
        unique_ptr<lmap_file, void (*)(lmap_file *)> f(lmap_open_io(store->openMember(g.hgnPath), 1), lmap_close);
        if (!f) { fail(gname + ": cannot open " + g.hgnPath); continue; }
        lmap_file *lf = f.get();
        if (strcmp(lmap_profile(lf), hal3::HGN_PROFILE) != 0) fail(gname + ": profile is " + lmap_profile(lf));
        const char *mg = lmap_meta_get(lf, "hal.genome"), *mp = lmap_meta_get(lf, "hal.parent");
        if (!mg || g.name != mg) fail(gname + ": the .hgn says it is " + (mg ? mg : "(none)"));
        if (g.parent != (mp ? mp : "")) fail(gname + ": the .hgn names parent " + (mp ? mp : "(none)"));
        if (lmap_verify(lf) != 0) fail(gname + ": liftmap integrity check failed");
        if (lmap_a_overlap(lf)) fail(gname + ": axis a overlaps -- some base has more than one parent position");
        lmap_file_stats st;
        lmap_get_stats(lf, &st);
        if (g.parent.empty() && (st.runs || lmap_n_members(lf, 1)))
            fail(gname + ": the root has runs or parent sequences");
        vector<pair<string, uint64_t>> &t = tables[g.name];
        uint64_t total = 0;
        for (uint32_t i = 0; i < lmap_n_members(lf, 0); i++) {
            const lmap_member *mm = lmap_member_at(lf, 0, i);
            t.push_back({mm->name, mm->length});
            total += mm->length;
        }
        if (t.size() != g.numSequences || total != g.length)
            fail(gname + ": " + to_string(t.size()) + " sequences, " + to_string(total) + " bp; the manifest says " +
                 to_string(g.numSequences) + ", " + to_string(g.length));
        if (!g.dnaPath.empty()) {
            unique_ptr<tb_file, void (*)(tb_file *)> d(nullptr, tb_close);
            if (!store->has(g.dnaPath)) fail(gname + ": " + g.dnaPath + " is missing");
            else {
                lmap_io *io = store->openMember(g.dnaPath);
                d.reset(tb_open_io(io, 1));
                if (!d) {
                    lmap_io_close(io);
                    fail(gname + ": cannot open " + g.dnaPath);
                }
            }
            if (d) {
                if (tb_n_seqs(d.get()) != t.size()) fail(gname + ": the .2bit has " + to_string(tb_n_seqs(d.get())) + " sequences");
                else
                    for (uint32_t i = 0; i < t.size(); i++) {
                        char err[256];
                        if (t[i].first != tb_name(d.get(), i) || tb_length(d.get(), i) != (int64_t)t[i].second) {
                            fail(gname + ": .2bit sequence " + to_string(i) + " is " + tb_name(d.get(), i) + ", not " + t[i].first);
                            break;
                        }
                        if (tb_verify_seq(d.get(), i, err, sizeof err) != 0) fail(gname + ": .2bit " + t[i].first + ": " + err);
                    }
            }
        }
        /* likelihoods: for this genome's sequences, every chunk sound */
        if (const hal3::MemberRef *mr = g.member("lik")) {
            if (store->has(mr->path)) {
                try {
                    unique_ptr<hal3::LikFile> lk = hal3::LikFile::open(store->openMember(mr->path), true, mr->path);
                    bool same = lk->numSequences() == t.size();
                    for (uint32_t i = 0; same && i < t.size(); i++)
                        same = lk->sequenceName(i) == t[i].first && lk->sequenceLength(i) == (int64_t)t[i].second;
                    if (!same) fail(gname + ": " + mr->path + " is not for its sequences");
                    else
                        for (const string &p : lk->verify()) fail(gname + ": " + p);
                    if (g.dnaPath.empty()) fail(gname + ": has likelihoods but no DNA (they are relative to its bases)");
                } catch (exception &e) {
                    fail(gname + ": " + e.what());
                }
            }
        }
        /* axis b against the parent's axis a (the parent comes first in the roster) */
        if (!g.parent.empty()) {
            auto pt = tables.find(g.parent);
            if (pt == tables.end()) fail(gname + ": listed before its parent " + g.parent);
            else {
                uint32_t nb = lmap_n_members(lf, 1);
                if (nb != pt->second.size()) fail(gname + ": axis b has " + to_string(nb) + " sequences, the parent " + to_string(pt->second.size()));
                else
                    for (uint32_t i = 0; i < nb; i++) {
                        const lmap_member *mm = lmap_member_at(lf, 1, i);
                        if (pt->second[i].first != mm->name || pt->second[i].second != mm->length) {
                            fail(gname + ": parent sequence " + mm->name + " differs from " + g.parent + "'s own table");
                            break;
                        }
                    }
            }
        }
    }

    /* ---- members the manifest does not name; in a pack, members against its index */
    for (const string &n : store->members())
        if (n != hal3::MANIFEST_NAME && !named.count(n)) fail(n + " is not in the manifest");
    if (store->isPack() && doCrc) {
        for (const auto &e : hal3::packInfo(*store).entries)
            if (store->crc(e.first) != e.second.crc) fail(e.first + ": CRC differs from the pack's index");
    }

    cout << dir << ": " << m.genomes.size() << " genomes, " << (errors ? "FAILED (" + to_string(errors) + " problems)" : "ok") << endl;
    return errors ? 1 : 0;
}

int main_validate(int argc, char **argv) { return runCommand("validate", [&]() { return validate(argc, argv); }); }
}
