/*
 * halToFormat3: convert a HAL alignment to HAL format 3 (the unpacked directory form).
 *
 * Genomes are independent -- each reads only itself and its parent -- so they convert
 * in parallel, in worker processes (HDF5 is not thread-safe), each opening the input
 * itself.  The manifest is written last, so an interrupted run leaves no manifest and
 * is not mistaken for a finished alignment.
 */
#include "hal.h"
#include "hal3Convert.h"
#include "hal3Manifest.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace std;
using namespace hal;

static void initParser(CLParser *optionsParser) {
    optionsParser->addArgument("inHalPath", "input HAL file");
    optionsParser->addArgument("outDir", "output directory (the format 3 alignment); must not already hold a manifest");
    optionsParser->addOption("threads", "genomes to convert at once (worker processes)", 1);
    optionsParser->addOption("check", "random bases per genome to check against the input (0: none)", 10000);
    optionsParser->addOption("chunkRuns", "at most this many runs per chunk", 8192);
    optionsParser->addOption("chunkParentSpan", "at most this much parent sequence per chunk", 1000000);
    optionsParser->addOptionFlag("noDna", "do not write .2bit files", false);
    optionsParser->addOptionFlag("checkOnly", "convert nothing: check an existing outDir against inHalPath", false);
    optionsParser->setDescription("Convert a HAL alignment to HAL format 3: a directory holding a manifest, and per "
                                  "genome a .hgn genome file (its branch to the parent) and a .2bit of its DNA.");
}

static void mkdirOrDie(const string &path) {
    if (mkdir(path.c_str(), 0777) != 0 && errno != EEXIST)
        throw hal_exception("cannot create " + path + ": " + strerror(errno));
}

/* Preorder from the root: the manifest's order, parents before children. */
static void preorder(const Alignment *aln, const string &name, vector<string> &out) {
    out.push_back(name);
    for (const string &c : aln->getChildNames(name)) preorder(aln, c, out);
}

int main(int argc, char **argv) {
    CLParser optionsParser;
    initParser(&optionsParser);
    string inPath, outDir;
    int threads;
    uint64_t checkSamples;
    bool checkOnly;
    hal3::ConvertOptions opts;
    try {
        optionsParser.parseOptions(argc, argv);
        inPath = optionsParser.getArgument<string>("inHalPath");
        outDir = optionsParser.getArgument<string>("outDir");
        threads = max(1, optionsParser.getOption<int>("threads"));
        checkSamples = optionsParser.getOption<uint64_t>("check");
        opts.chunkRuns = optionsParser.getOption<uint32_t>("chunkRuns");
        opts.chunkParentSpan = optionsParser.getOption<uint64_t>("chunkParentSpan");
        opts.dna = !optionsParser.getFlag("noDna");
        checkOnly = optionsParser.getFlag("checkOnly");
    } catch (exception &e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        return 1;
    }

    try {
        string manifestPath = outDir + "/" + hal3::MANIFEST_NAME;
        if (!checkOnly && access(manifestPath.c_str(), F_OK) == 0)
            throw hal_exception(outDir + " already holds an alignment (" + manifestPath + ")");
        if (checkOnly && !checkSamples) checkSamples = 10000;
        string genomesDir = outDir + "/genomes";

        /* the tree, roster and metadata; the input is closed again before forking */
        hal3::Manifest manifest;
        vector<string> names;
        {
            AlignmentConstPtr aln(openHalAlignment(inPath, &optionsParser));
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
                AlignmentConstPtr aln(openHalAlignment(inPath, &optionsParser));
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
                cerr << "halToFormat3: " << e.what() << endl;
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
    } catch (exception &e) {
        cerr << "halToFormat3: " << e.what() << endl;
        return 1;
    }
    return 0;
}
