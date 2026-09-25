/*
 * hal coverage: how much of a reference genome each leaf genome covers, and how many
 * times, as halCoverage.  A reference base is covered k times by a genome when liftover
 * to that genome -- paralogs included, through their common ancestor, as halLiftover and
 * halCoverage map -- gives k bases.  Exact over every base by default; --numSamples
 * samples positions as halCoverage does.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Walk.h"

#include <algorithm>
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <map>
#include <random>

using namespace std;

namespace hal3 {

static const char *COVERAGE_USAGE =
    "usage: hal coverage [options] <alignment> <refGenome>\n"
    "Coverage of refGenome by each leaf genome, as halCoverage: for each leaf, the reference\n"
    "bases it covers at least 1, 2, ... times (liftover, paralogs included).  Every base is\n"
    "counted unless --numSamples is given.\n"
    "options:\n"
    "  --numSamples N   count N random reference positions, drawn with replacement, as\n"
    "                   halCoverage does, instead of every base\n"
    "  --seed S         random seed for --numSamples (default: the time)\n"
    "  --bySequence     also a table for each reference sequence\n"
    "  --refSequence S  only this reference sequence (and so only its table)\n"
    "  --threads N      count N chunks of the reference at once (default 1)\n"
    "  -h, --help\n";

namespace {

/* Cumulative histogram: h[k] = positions covered at least k+1 times. */
typedef vector<uint64_t> Hist;

void add(Hist &h, int depth, uint64_t weight) {
    if (depth <= 0) return;
    if ((int)h.size() < depth) h.resize((size_t)depth, 0);
    for (int k = 0; k < depth; k++) h[(size_t)k] += weight;
}

void add(Hist &h, const Hist &o) {
    if (h.size() < o.size()) h.resize(o.size(), 0);
    for (size_t k = 0; k < o.size(); k++) h[k] += o[k];
}

/* The depth of pieces over [cs, ce), as (position, depth) steps in order. */
vector<pair<int64_t, int>> depthSteps(const vector<LiftPiece> &pieces) {
    vector<pair<int64_t, int>> ev;
    for (const LiftPiece &p : pieces) {
        ev.push_back({p.srcStart, 1});
        ev.push_back({p.srcStart + p.length, -1});
    }
    sort(ev.begin(), ev.end());
    vector<pair<int64_t, int>> steps;
    int d = 0;
    for (size_t i = 0; i < ev.size();) {
        int64_t x = ev[i].first;
        for (; i < ev.size() && ev[i].first == x; i++) d += ev[i].second;
        steps.push_back({x, d});
    }
    return steps;
}
}

int main_coverage(int argc, char **argv) {
    return runCommand("coverage", [&]() {
        int64_t numSamples = 0, seed = 0, threads = 1, chunk = 2000000;
        bool bySequence = false;
        string refSeqName;
        enum { NUM_SAMPLES = 1000, SEED, BY_SEQUENCE, REF_SEQUENCE, THREADS, CHUNK };
        static const option longOpts[] = {{"numSamples", required_argument, nullptr, NUM_SAMPLES},
                                          {"seed", required_argument, nullptr, SEED},
                                          {"bySequence", no_argument, nullptr, BY_SEQUENCE},
                                          {"refSequence", required_argument, nullptr, REF_SEQUENCE},
                                          {"threads", required_argument, nullptr, THREADS},
                                          {"chunk", required_argument, nullptr, CHUNK},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case NUM_SAMPLES: numSamples = parseInt("numSamples", optarg); break;
            case SEED: seed = parseInt("seed", optarg); break;
            case BY_SEQUENCE: bySequence = true; break;
            case REF_SEQUENCE: refSeqName = optarg; break;
            case THREADS: threads = parseInt("threads", optarg); break;
            case CHUNK: chunk = parseInt("chunk", optarg); break;
            case 'h': cout << COVERAGE_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        if (numSamples < 0) throw CliError("--numSamples must not be negative");
        if (threads < 1 || chunk < 1) throw CliError("--threads and --chunk must be positive");
        Alignment aln(argv[optind]);
        const Genome *ref = aln.genome(argv[optind + 1]);
        if (!ref) throw runtime_error(string("Genome ") + argv[optind + 1] + " not found");

        /* leaves in tree order */
        vector<const Genome *> leaves, stack{aln.root()};
        while (!stack.empty()) {
            const Genome *g = stack.back();
            stack.pop_back();
            if (g->children().empty()) leaves.push_back(g);
            for (auto ch = g->children().rbegin(); ch != g->children().rend(); ++ch) stack.push_back(*ch);
        }
        WalkOptions wo;
        wo.onlyOrthologs = true;
        for (const Genome *l : leaves) wo.targets.push_back(l->name());
        { TreeWalk check(aln, ref, wo); }

        vector<uint32_t> seqs;
        if (!refSeqName.empty()) {
            int64_t i = ref->findSequence(refSeqName);
            if (i < 0) throw runtime_error("no sequence " + refSeqName + " in " + ref->name());
            seqs.push_back((uint32_t)i);
        } else {
            for (uint32_t i = 0; i < ref->numSequences(); i++) seqs.push_back(i);
        }
        /* sampled positions, by sequence, sorted, with their multiplicity */
        map<uint32_t, map<int64_t, uint64_t>> samples;
        if (numSamples > 0) {
            int64_t total = 0;
            for (uint32_t i : seqs) total += ref->sequenceLength(i);
            if (total == 0) throw runtime_error("nothing to sample");
            mt19937_64 rng(seed ? (uint64_t)seed : (uint64_t)time(nullptr));
            uniform_int_distribution<int64_t> pick(0, total - 1);
            for (int64_t k = 0; k < numSamples; k++) {
                int64_t x = pick(rng);
                for (uint32_t i : seqs) {
                    if (x < ref->sequenceLength(i)) {
                        samples[i][x]++;
                        break;
                    }
                    x -= ref->sequenceLength(i);
                }
            }
        }

        /* chunks, counted on their own (each with its own walk, so its own cursors) and
         * added up in order */
        vector<RefInterval> whole0;
        for (uint32_t si : seqs) whole0.push_back(RefInterval{si, 0, ref->sequenceLength(si)});
        const vector<RefInterval> chunks = chunkIntervals(whole0, chunk);
        map<const Genome *, Hist> whole;
        map<uint32_t, map<const Genome *, Hist>> bySeq;
        orderedParallel<map<const Genome *, Hist>>(
            chunks.size(), (int)threads,
            [&](size_t k) {
                const uint32_t si = chunks[k].seq;
                const int64_t cs = chunks[k].start, ce = chunks[k].end;
                map<const Genome *, Hist> part;
                vector<pair<int64_t, int64_t>> intervals;
                vector<pair<int64_t, uint64_t>> points;   /* sampled positions in this chunk */
                if (numSamples > 0) {
                    auto sm = samples.find(si);
                    if (sm != samples.end())
                        for (auto it = sm->second.lower_bound(cs); it != sm->second.end() && it->first < ce; ++it) {
                            points.push_back(*it);
                            if (!intervals.empty() && intervals.back().second == it->first) intervals.back().second++;
                            else intervals.push_back({it->first, it->first + 1});
                        }
                    if (intervals.empty()) return part;
                } else {
                    intervals.push_back({cs, ce});
                }
                TreeWalk walk(aln, ref, wo);
                map<const Genome *, vector<LiftPiece>> pieces;
                walk.lift(si, intervals, pieces);
                for (const Genome *g : leaves) {
                    auto pg = pieces.find(g);
                    vector<pair<int64_t, int>> steps = depthSteps(pg == pieces.end() ? vector<LiftPiece>() : pg->second);
                    /* the reference covers itself once, as liftover to itself gives the base */
                    int self = g == ref ? 1 : 0;
                    Hist &h = part[g];
                    if (numSamples > 0) {
                        size_t i = 0;
                        int d = 0;
                        for (auto &pt : points) {
                            while (i < steps.size() && steps[i].first <= pt.first) d = steps[i++].second;
                            add(h, d + self, pt.second);
                        }
                    } else {
                        int64_t x = cs;
                        int d = 0;
                        for (auto &st : steps) {
                            if (st.first > x) {
                                add(h, d + self, (uint64_t)(st.first - x));
                                x = st.first;
                            }
                            d = st.second;
                        }
                        if (ce > x) add(h, d + self, (uint64_t)(ce - x));
                    }
                }
                return part;
            },
            [&](size_t k, map<const Genome *, Hist> &part) {
                for (auto &gh : part) {
                    add(whole[gh.first], gh.second);
                    if (bySequence) add(bySeq[chunks[k].seq][gh.first], gh.second);
                }
            });

        size_t maxDepth = 0;
        for (auto &w : whole) maxDepth = max(maxDepth, w.second.size());
        cout << "Genome";
        for (size_t i = 0; i < maxDepth; i++) cout << ", sitesCovered" << i + 1 << "Times";
        cout << endl;
        auto table = [&](map<const Genome *, Hist> &t) {
            for (const Genome *g : leaves) {
                const Hist &h = t[g];
                cout << g->name();
                for (size_t i = 0; i < maxDepth; i++) cout << ", " << (i < h.size() ? h[i] : 0);
                cout << endl;
            }
        };
        if (refSeqName.empty()) table(whole);
        if (bySequence || !refSeqName.empty()) {
            for (uint32_t si : seqs) {
                if (bySequence || refSeqName.empty()) cout << "\nCoverage on " << ref->sequenceName(si) << endl;
                table(bySequence ? bySeq[si] : whole);
            }
        }
        cout.flush();
        if (!cout) throw runtime_error("write failed");
        return 0;
    });
}
}
