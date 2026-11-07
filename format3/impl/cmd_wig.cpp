/*
 * hal wig: scores on the universal axis (hal maf --universal --axisRow), projected onto a
 * genome -- each base gets its column's score.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Liftover.h"
#include "hal3Universal.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>

using namespace std;

namespace hal3 {

static const char *WIG_USAGE =
    "usage: hal wig [options] <alignment> <axis.wig> <genome> <out.wig>\n"
    "Scores of universal columns -- a wig on the axis of hal maf --universal --axisRow, as\n"
    "phyloP --wig-scores gives for it -- projected onto a genome: each of its bases gets the\n"
    "score of the column it is in.  The axis options must be the export's.  A base in no\n"
    "column (a leaf's base with no parent, without --leafInsertions; with --noDupes, a\n"
    "leaf's copy that is not the canonical one) or in a column with no score gets none.\n"
    "out.wig may be stdout; its sequences are named as the genome's (--ucscSequenceNames:\n"
    "genome.sequence).  Memory: 4 bytes per column of the axis sequences the wig names, and\n"
    "per base of the genome's longest sequence.\n"
    "options (as the export's):\n"
    "  --noDupes            columns with one copy per genome\n"
    "  --leafInsertions     leaves' bases with no parent have columns\n"
    "  --rootGenome R       the axis's backbone is R, not the root\n"
    "  --excludeGenomes a,b these leaves have no rows\n"
    "options:\n"
    "  --decimals N         digits after the point (default 3, as phyloP)\n"
    "  --ucscSequenceNames  name sequences genome.sequence\n"
    "  --threads N          threads building the axis (default 1)\n"
    "  -h, --help\n";

namespace {

/* A run of the axis: bases [pos, pos+len) of a sequence of a genome at columns
 * [col, col+len) of an axis sequence (decreasing bases when rev). */
struct Run {
    int64_t code;   /* of its first base, pos */
    int64_t pos, len, col;
    uint32_t range;
    bool rev;
};

/* The wig's scores, by axis sequence and column; NaN where it has none.  An axis
 * sequence the wig does not name is left empty (no memory). */
void readWig(const string &path, const map<string, uint32_t> &rangeOf, const vector<int64_t> &rangeCols,
             vector<vector<float>> &scores) {
    ifstream in(path == "-" ? "/dev/stdin" : path);
    if (!in) throw runtime_error("cannot open " + path);
    string line;
    int64_t pos = -1, step = 1, span = 1;
    bool fixed = true;
    vector<float> *cur = nullptr;
    set<string> unknown;
    uint64_t lineNo = 0;
    while (getline(in, line)) {
        lineNo++;
        if (line.empty() || line[0] == '#' || !line.compare(0, 5, "track") || !line.compare(0, 7, "browser")) continue;
        if (!line.compare(0, 9, "fixedStep") || !line.compare(0, 12, "variableStep")) {
            fixed = line[0] == 'f';
            map<string, string> kv;
            stringstream ss(line);
            string w;
            ss >> w;
            while (ss >> w) {
                size_t eq = w.find('=');
                if (eq != string::npos) kv[w.substr(0, eq)] = w.substr(eq + 1);
            }
            auto r = rangeOf.find(kv["chrom"]);
            if (r == rangeOf.end()) {
                unknown.insert(kv["chrom"]);
                cur = nullptr;
            } else {
                cur = &scores[r->second];
                if (cur->empty()) cur->assign((size_t)rangeCols[r->second], numeric_limits<float>::quiet_NaN());
            }
            pos = fixed ? stoll(kv.count("start") ? kv["start"] : "1") - 1 : -1;
            step = kv.count("step") ? stoll(kv["step"]) : 1;
            span = kv.count("span") ? stoll(kv["span"]) : 1;
            continue;
        }
        if (!cur && unknown.empty()) throw runtime_error(path + " line " + to_string(lineNo) + ": data before a fixedStep or variableStep line");
        int64_t at;
        double v;
        if (fixed) {
            v = atof(line.c_str());
            at = pos;
            pos += step;
        } else {
            char *end;
            at = strtoll(line.c_str(), &end, 10) - 1;
            v = atof(end);
        }
        if (!cur) continue;
        for (int64_t k = at; k < at + span; k++) {
            if (k < 0 || k >= (int64_t)cur->size())
                throw runtime_error(path + " line " + to_string(lineNo) + ": position " + to_string(k + 1) + " is past its axis sequence");
            (*cur)[(size_t)k] = (float)v;
        }
    }
    if (!unknown.empty())
        throw runtime_error(path + ": " + *unknown.begin() + (unknown.size() > 1 ? " (and others)" : "") +
                            " is not an axis sequence of this alignment with these options");
}
}

int main_wig(int argc, char **argv) {
    return runCommand("wig", [&]() {
        bool dupes = true, leafInsertions = false, ucsc = false;
        string rootGenome;
        vector<string> exclude;
        int64_t decimals = 3, threads = 1;
        enum { NO_DUPES = 1000, LEAF_INSERTIONS, ROOT_GENOME, EXCLUDE_GENOMES, DECIMALS, UCSC, THREADS };
        static const option longOpts[] = {{"noDupes", no_argument, nullptr, NO_DUPES},
                                          {"leafInsertions", no_argument, nullptr, LEAF_INSERTIONS},
                                          {"rootGenome", required_argument, nullptr, ROOT_GENOME},
                                          {"excludeGenomes", required_argument, nullptr, EXCLUDE_GENOMES},
                                          {"decimals", required_argument, nullptr, DECIMALS},
                                          {"ucscSequenceNames", no_argument, nullptr, UCSC},
                                          {"threads", required_argument, nullptr, THREADS},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        optind = 1;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case NO_DUPES: dupes = false; break;
            case LEAF_INSERTIONS: leafInsertions = true; break;
            case ROOT_GENOME: rootGenome = optarg; break;
            case EXCLUDE_GENOMES: {
                stringstream ss(optarg);
                string x;
                while (getline(ss, x, ','))
                    if (!x.empty()) exclude.push_back(x);
                break;
            }
            case DECIMALS: decimals = parseInt("decimals", optarg); break;
            case UCSC: ucsc = true; break;
            case THREADS: threads = parseInt("threads", optarg); break;
            case 'h': cout << WIG_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 4) throw CliError("expected 4 arguments, got " + to_string(argc - optind));
        if (decimals < 0 || decimals > 9 || threads < 1) throw CliError("--decimals must be 0 to 9, --threads positive");
        Alignment aln(argv[optind]);
        const string wigPath = argv[optind + 1], outPath = argv[optind + 3];
        const Genome *target = aln.genome(argv[optind + 2]);
        if (!target) throw runtime_error(string("no genome ") + argv[optind + 2]);
        const Genome *backbone = rootGenome.empty() ? aln.root() : aln.genome(rootGenome);
        if (!backbone) throw runtime_error("no genome " + rootGenome);
        bool below = false;
        for (const Genome *a = target; a; a = a->parent()) below = below || a == backbone;
        if (!below) throw runtime_error(target->name() + " is not in the subtree of " + backbone->name());
        set<const Genome *> excluded;
        for (const string &x : exclude) {
            const Genome *g = aln.genome(x);
            if (!g) throw runtime_error("no genome " + x);
            if (!g->children().empty()) throw CliError("only leaves can be excluded (" + x + " is an ancestor)");
            excluded.insert(g);
        }
        if (excluded.count(target)) throw CliError(target->name() + " is excluded: it has no rows");

        /* the axis, every column of it: the backbone's sequences, then the flankless */
        UniversalIndex index(aln, backbone, leafInsertions, dupes, excluded, false, (int)threads);
        vector<TopRange> ranges;
        for (uint32_t q = 0; q < backbone->numSequences(); q++) ranges.push_back(TopRange{0, q, 0, backbone->sequenceLength(q)});
        UniversalPlanner plan(index, ranges, true, (int64_t)1 << 30);
        plan.run();
        const size_t nRanges = plan.ranges().size();
        vector<int64_t> rangeU0(nRanges, -1), rangeCols(nRanges, 0);
        for (size_t k = 0; k < plan.chunks(); k++) {
            const UniversalPlanner::Chunk &ch = plan.chunk(k);
            if (rangeU0[ch.range] < 0) rangeU0[ch.range] = ch.u0;
            rangeCols[ch.range] += ch.cols;
        }
        map<string, uint32_t> rangeOf;
        for (size_t r = 0; r < nRanges; r++) {
            const TopRange &t = plan.ranges()[r];
            const Genome *g = index.genome(t.g);
            rangeOf["_axis." + g->name() + "." + g->sequenceName(t.seq)] = (uint32_t)r;
        }
        /* the runs of top bases, by code, of the genomes a base of the target can top out
         * in: it and its ancestors */
        vector<vector<Run>> runs((size_t)index.genomes());
        vector<char> onPath((size_t)index.genomes(), 0);
        for (const Genome *a = target; a; a = a == backbone ? nullptr : a->parent()) onPath[(size_t)index.rank(a)] = 1;
        vector<Segment> segs;
        for (size_t k = 0; k < plan.chunks(); k++) {
            const UniversalPlanner::Chunk &ch = plan.chunk(k);
            plan.segments(k, segs);
            for (const Segment &sg : segs)
                if (onPath[(size_t)sg.g]) runs[(size_t)sg.g].push_back(Run{index.code(sg.g, sg.seq, sg.pos), sg.pos, sg.len, ch.u0 - rangeU0[ch.range] + sg.u,
                                                 (uint32_t)ch.range, sg.rev});
        }
        for (auto &v : runs) sort(v.begin(), v.end(), [](const Run &a, const Run &b) { return a.code < b.code; });

        vector<vector<float>> scores(nRanges);
        readWig(wigPath, rangeOf, rangeCols, scores);

        FILE *out = outPath == "stdout" || outPath == "-" ? stdout : fopen(outPath.c_str(), "w");
        if (!out) throw runtime_error("cannot write " + outPath);
        CursorCache cursors;
        vector<float> seqScores;
        vector<LiftPiece> cur, next, unmapped;
        vector<pair<const Genome *, vector<LiftPiece>>> tops;
        uint64_t projected = 0, inColumns = 0;
        char buf[64];
        for (uint32_t q = 0; q < target->numSequences(); q++) {
            const int64_t len = target->sequenceLength(q);
            seqScores.assign((size_t)len, numeric_limits<float>::quiet_NaN());
            /* up to each base's top: where its column comes from */
            tops.clear();
            cur.assign(1, LiftPiece{0, q, 0, len, false});
            for (const Genome *a = target; !cur.empty(); a = a->parent()) {
                if (a == backbone) {
                    tops.emplace_back(a, cur);
                    break;
                }
                unmapped.clear();
                liftUpBranch(a, cur, next, &unmapped, cursors, !dupes);
                if (!unmapped.empty()) tops.emplace_back(a, unmapped);
                cur.swap(next);
            }
            for (auto &t : tops) {
                const int rank = index.rank(t.first);
                if (rank < 0 || (rank > 0 && !index.covered(rank))) continue;   /* the backbone: every base a column */
                const vector<Run> &rv = runs[(size_t)rank];
                for (const LiftPiece &p : t.second) {
                    /* the runs overlapping the top bases [tgtStart, +length) of sequence tgtSeq */
                    const int64_t lo = index.code(rank, p.tgtSeq, p.tgtStart), hi = lo + p.length;
                    auto it = upper_bound(rv.begin(), rv.end(), lo, [](int64_t x, const Run &r) { return x < r.code; });
                    if (it != rv.begin()) --it;
                    for (; it != rv.end() && it->code < hi; ++it) {
                        const int64_t a = max(lo, it->code), b = min(hi, it->code + it->len);
                        if (a >= b) continue;
                        const vector<float> &sc = scores[it->range];
                        if (sc.empty()) continue;   /* an axis sequence the wig has no score on */
                        for (int64_t x = a; x < b; x++) {
                            const int64_t i = x - lo;   /* the piece's i-th top base */
                            const int64_t pos = p.reversed ? p.srcStart + p.length - 1 - i : p.srcStart + i;
                            const int64_t k = x - it->code;   /* the run's k-th base */
                            const int64_t col = it->rev ? it->col + it->len - 1 - k : it->col + k;
                            seqScores[(size_t)pos] = sc[(size_t)col];
                            inColumns++;
                        }
                    }
                }
            }
            /* the sequence's scores, as fixedStep runs */
            const string name = ucsc ? target->name() + "." + target->sequenceName(q) : string(target->sequenceName(q));
            bool open = false;
            for (int64_t x = 0; x < len; x++) {
                const float v = seqScores[(size_t)x];
                if (std::isnan(v)) {
                    open = false;
                    continue;
                }
                if (!open) {
                    fprintf(out, "fixedStep chrom=%s start=%lld step=1\n", name.c_str(), (long long)(x + 1));
                    open = true;
                }
                snprintf(buf, sizeof buf, "%.*f\n", (int)decimals, (double)v);
                fputs(buf, out);
                projected++;
            }
        }
        if (out != stdout && fclose(out) != 0) throw runtime_error("cannot write " + outPath);
        cerr << "hal wig: " << target->name() << ": " << inColumns << " bases in columns of the wig's axis sequences, " << projected << " with a score"
             << endl;
        return 0;
    });
}
}
