/*
 * hal convert: a HAL (HDF5 or mmap) alignment to HAL format 3 (the unpacked directory
 * form).  The one command that needs HDF5, and so libHal.
 *
 * Genomes are independent -- each reads only itself and its parent -- so they convert
 * in parallel, in worker processes (HDF5 is not thread-safe), each opening the input
 * itself.  The manifest is written last, so an interrupted run leaves no manifest and
 * is not mistaken for a finished alignment.
 */
#include "hal.h"
#include "hal3Cli.h"
#include "hal3Convert.h"
#include "hal3Manifest.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace std;
using namespace hal;

static const char *CONVERT_USAGE =
    "usage: hal convert [options] <in.hal> <outDir>\n"
    "Convert a HAL (HDF5 or mmap) alignment to format 3: a directory holding a manifest, and per\n"
    "genome a .hgn genome file (its branch to the parent) and a .2bit of its DNA.  outDir must\n"
    "not already hold an alignment.  Every genome is checked against the input afterwards.\n"
    "options:\n"
    "  --threads N          genomes to convert at once, in worker processes (default 1)\n"
    "  --check N            random bases per genome to check against the input (default 10000;\n"
    "                       0: none)\n"
    "  --checkOnly          convert nothing: check an existing outDir against in.hal\n"
    "  --chunkRuns N        at most N runs per chunk (default 8192)\n"
    "  --chunkParentSpan N  at most N bases of parent sequence per chunk (default 1000000)\n"
    "  --noDna              do not write .2bit files\n"
    "  -h, --help\n";

static void mkdirOrDie(const string &path) {
    if (mkdir(path.c_str(), 0777) != 0 && errno != EEXIST)
        throw hal_exception("cannot create " + path + ": " + strerror(errno));
}

/* Preorder from the root: the manifest's order, parents before children. */
static void preorder(const Alignment *aln, const string &name, vector<string> &out) {
    out.push_back(name);
    for (const string &c : aln->getChildNames(name)) preorder(aln, c, out);
}

static int convert(int argc, char **argv) {
    string inPath, outDir;
    int threads = 1;
    uint64_t checkSamples = 10000;
    bool checkOnly = false;
    hal3::ConvertOptions opts;
    {
        enum { THREADS = 1000, CHECK, CHECK_ONLY, CHUNK_RUNS, CHUNK_PARENT_SPAN, NO_DNA };
        static const option longOpts[] = {{"threads", required_argument, nullptr, THREADS},
                                          {"check", required_argument, nullptr, CHECK},
                                          {"checkOnly", no_argument, nullptr, CHECK_ONLY},
                                          {"chunkRuns", required_argument, nullptr, CHUNK_RUNS},
                                          {"chunkParentSpan", required_argument, nullptr, CHUNK_PARENT_SPAN},
                                          {"noDna", no_argument, nullptr, NO_DNA},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case THREADS: threads = (int)max<int64_t>(1, hal3::parseInt("threads", optarg)); break;
            case CHECK: checkSamples = (uint64_t)max<int64_t>(0, hal3::parseInt("check", optarg)); break;
            case CHECK_ONLY: checkOnly = true; break;
            case CHUNK_RUNS: opts.chunkRuns = (uint32_t)hal3::parseInt("chunkRuns", optarg); break;
            case CHUNK_PARENT_SPAN: opts.chunkParentSpan = (uint64_t)hal3::parseInt("chunkParentSpan", optarg); break;
            case NO_DNA: opts.dna = false; break;
            case 'h': cout << CONVERT_USAGE; return 0;
            default: throw hal3::CliError("bad option");
            }
        }
        if (argc - optind != 2) throw hal3::CliError("expected 2 arguments, got " + to_string(argc - optind));
        inPath = argv[optind];
        outDir = argv[optind + 1];
    }

    {
        string manifestPath = outDir + "/" + hal3::MANIFEST_NAME;
        if (!checkOnly && access(manifestPath.c_str(), F_OK) == 0)
            throw hal_exception(outDir + " already holds an alignment (" + manifestPath + ")");
        if (checkOnly && !checkSamples) checkSamples = 10000;
        string genomesDir = outDir + "/genomes";

        /* the tree, roster and metadata; the input is closed again before forking */
        hal3::Manifest manifest;
        vector<string> names;
        {
            AlignmentConstPtr aln(openHalAlignment(inPath, nullptr));
            manifest.newick = aln->getNewickTree();
            try {
                for (const auto &kv : aln->getMetaData()->getMap()) manifest.meta[kv.first] = kv.second;
            } catch (hal_exception &) {
                /* the mmap backend keeps no alignment metadata */
            }
            preorder(aln.get(), aln->getRootName(), names);
            for (const string &n : names) {
                const Genome *g = aln->openGenomeCheck(n);
                hal3::GenomeEntry e;
                e.name = n;
                e.parent = aln->getParentName(n);
                e.branchLength = e.parent.empty() ? -1 : aln->getBranchLength(e.parent, n);
                e.length = g->getSequenceLength();
                e.numSequences = g->getNumSequences();
                e.hgnPath = "genomes/" + hal3::genomeFileStem(n) + ".hgn";
                if (opts.dna && g->containsDNAArray()) e.dnaPath = "genomes/" + hal3::genomeFileStem(n) + ".2bit";
                manifest.genomes.push_back(e);
                aln->closeGenome(g);
            }
        }
        if (!checkOnly) {
            mkdirOrDie(outDir);
            mkdirOrDie(genomesDir);
        }

        /* largest genomes first, dealt round-robin to the workers */
        vector<size_t> order(names.size());
        for (size_t i = 0; i < order.size(); i++) order[i] = i;
        stable_sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return manifest.genomes[a].length > manifest.genomes[b].length; });
        int nw = (int)min<size_t>((size_t)threads, names.size());
        vector<pid_t> pids;
        for (int w = 0; w < nw; w++) {
            cout.flush();
            cerr.flush();
            pid_t pid = fork();
            if (pid < 0) throw hal_exception(string("fork failed: ") + strerror(errno));
            if (pid > 0) {
                pids.push_back(pid);
                continue;
            }
            int rc = 0;
            try {
                AlignmentConstPtr aln(openHalAlignment(inPath, nullptr));
                for (size_t k = (size_t)w; k < order.size(); k += (size_t)nw) {
                    const string &n = names[order[k]];
                    hal3::GenomeConversion c;
                    if (!checkOnly) c = hal3::convertGenome(aln.get(), n, genomesDir, opts);
                    uint64_t bad = checkSamples ? hal3::checkGenome(aln.get(), n, genomesDir, checkSamples, k + 1, cerr) : 0;
                    if (checkOnly) {
                        cout << n << "\t" << (bad ? "CHECK FAILED (" + to_string(bad) + ")" : "checked") << endl;
                        if (bad) rc = 1;
                        continue;
                    }
                    cout << n << "\t" << c.topSegments << " top segments, " << c.alignedSegments << " aligned -> "
                         << c.runs << " runs"
                         << (c.alignedSegments ? " (" + to_string((double)c.alignedSegments / max<uint64_t>(1, c.runs)).substr(0, 4) + "x)" : "")
                         << (c.wroteDna ? ", DNA" : "")
                         << (checkSamples ? (bad ? ", CHECK FAILED (" + to_string(bad) + ")" : ", checked") : "") << endl;
                    if (bad) rc = 1;
                }
            } catch (exception &e) {
                cerr << "hal convert: " << e.what() << endl;
                rc = 1;
            }
            cout.flush();
            _exit(rc);
        }
        bool failed = false;
        for (pid_t p : pids) {
            int status = 0;
            if (waitpid(p, &status, 0) < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) failed = true;
        }
        if (failed) throw hal_exception(checkOnly ? "check failed" : "conversion failed; no manifest written");
        if (checkOnly) return 0;

        for (hal3::GenomeEntry &e : manifest.genomes) {
            e.hgnCrc = hal3::fileCrc32(outDir + "/" + e.hgnPath, &e.hgnBytes);
            if (!e.dnaPath.empty()) e.dnaCrc = hal3::fileCrc32(outDir + "/" + e.dnaPath, &e.dnaBytes);
        }
        manifest.write(manifestPath);
        cout << "wrote " << manifestPath << " (" << manifest.genomes.size() << " genomes)" << endl;
    }
    return 0;
}

int hal3::main_convert(int argc, char **argv) {
    return hal3::runCommand("convert", [&]() { return convert(argc, argv); });
}
