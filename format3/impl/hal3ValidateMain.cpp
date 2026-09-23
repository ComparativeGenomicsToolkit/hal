/*
 * hal3Validate: check a HAL format 3 alignment (the unpacked directory form).
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
 *   orphans     no .hgn or .2bit in genomes/ that the manifest does not name
 */
#include "hal3Manifest.h"
#include "liftmap.h"
#include "twobit64.h"

#include <cmath>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace std;

static int errors = 0;
static void fail(const string &what) {
    if (errors++ < 50) cerr << "  FAIL " << what << endl;
}

static void usage() {
    cerr << "usage: hal3Validate alignment_dir [--noCrc]\n"
            "Check a HAL format 3 alignment: manifest, tree, files, checksums, and every\n"
            "genome file against its parent's and its DNA.\n";
    exit(2);
}

int main(int argc, char **argv) {
    string dir;
    bool doCrc = true;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--noCrc")) doCrc = false;
        else if (argv[i][0] == '-') usage();
        else if (dir.empty()) dir = argv[i];
        else usage();
    }
    if (dir.empty()) usage();

    hal3::Manifest m;
    try {
        m = hal3::Manifest::read(dir + "/" + hal3::MANIFEST_NAME);
    } catch (exception &e) {
        cerr << "hal3Validate: " << e.what() << endl;
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
        if (doCrc) {
            try {
                uint64_t bytes;
                uint32_t crc = hal3::fileCrc32(dir + "/" + g.hgnPath, &bytes);
                if (bytes != g.hgnBytes || crc != g.hgnCrc) fail(gname + ": " + g.hgnPath + " size or CRC differs from the manifest");
                if (!g.dnaPath.empty()) {
                    crc = hal3::fileCrc32(dir + "/" + g.dnaPath, &bytes);
                    if (bytes != g.dnaBytes || crc != g.dnaCrc) fail(gname + ": " + g.dnaPath + " size or CRC differs from the manifest");
                }
            } catch (exception &e) {
                fail(gname + ": " + e.what());
            }
        }
        unique_ptr<lmap_file, void (*)(lmap_file *)> f(lmap_open((dir + "/" + g.hgnPath).c_str()), lmap_close);
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
            unique_ptr<tb_file, void (*)(tb_file *)> d(tb_open((dir + "/" + g.dnaPath).c_str()), tb_close);
            if (!d) fail(gname + ": cannot open " + g.dnaPath);
            else if (tb_n_seqs(d.get()) != t.size()) fail(gname + ": the .2bit has " + to_string(tb_n_seqs(d.get())) + " sequences");
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

    /* ---- orphaned files */
    if (DIR *d = opendir((dir + "/genomes").c_str())) {
        while (struct dirent *e = readdir(d)) {
            string n = e->d_name;
            bool ours = (n.size() > 4 && n.compare(n.size() - 4, 4, ".hgn") == 0) ||
                        (n.size() > 5 && n.compare(n.size() - 5, 5, ".2bit") == 0);
            if (ours && !named.count("genomes/" + n)) fail("genomes/" + n + " is not in the manifest");
        }
        closedir(d);
    }

    cout << dir << ": " << m.genomes.size() << " genomes, " << (errors ? "FAILED (" + to_string(errors) + " problems)" : "ok") << endl;
    return errors ? 1 : 0;
}
