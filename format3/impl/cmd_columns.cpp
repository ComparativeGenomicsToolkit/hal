/*
 * hal columns: universal columns as trees of bases (hal3Columns.h), as text -- to check
 * them, and to look at them.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Columns.h"

#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <unordered_map>

using namespace std;

namespace hal3 {

static const char *COLUMNS_USAGE =
    "usage: hal columns [options] <alignment> <out.txt>\n"
    "The universal columns (hal maf --universal, duplicates kept, every ancestor base) as\n"
    "trees of bases: each base with its parent, the base of the parent genome it aligns to,\n"
    "which is always in the same column.  One line per column:\n"
    "  axis-sequence <tab> position (0-based) <tab> entry <tab> entry ...\n"
    "with entries in genome preorder (a parent before its children), each\n"
    "  genome.sequence|position|strand|base|parent|piece\n"
    "position 0-based on the forward strand; base in the axis's orientation (complemented on\n"
    "the - strand); parent the index of its parent's entry in the line (0 is the first\n"
    "entry), - for the column's top; piece a number for the stretch of the genome the base\n"
    "lies in (two copies under one parent, with the same pieces over several columns, are\n"
    "one duplication).\n"
    "options:\n"
    "  --refSequence S      only this backbone sequence\n"
    "  --start N            only from this position of it (0-based; with --refSequence)\n"
    "  --length N           only this many bases (with --refSequence)\n"
    "  --rootGenome R       the subtree of R, with R as the backbone\n"
    "  --excludeGenomes A,B leave out these leaves\n"
    "  --leafInsertions     a leaf's bases with no parent are columns of their own\n"
    "  --noAncestorBases    ancestors' bases are N (their DNA is not read)\n"
    "  --chunk N            columns a chunk (default 1048576)\n"
    "  --threads N          chunks on N threads (default 1); the output is the same\n"
    "sites for phyloFit (a neutral model on any columns, e.g. a random sample):\n"
    "  --ss                 write a phast SS file (tuple size 1, a row per leaf, phast names)\n"
    "                       instead: each column cut where a base has two copies of a genome\n"
    "                       below it, into single-copy parts, each part with two leaves or more\n"
    "                       a site (as hal 4d --dupes split)\n"
    "  --sampleChunks K     only K chunks, evenly spaced over the axis\n"
    "  --sample N           keep N sites: those whose parts' roots hash lowest (as hal 4d)\n";

namespace {
uint64_t fnvHash(const string &s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return h;
}

char tupleBase(char c) {
    c = (char)toupper((unsigned char)c);
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' ? c : 'N';
}
}

int main_columns(int argc, char **argv) {
    return runCommand("columns", [&]() {
        ColumnOptions opts;
        string refSeqName;
        int64_t start = 0, length = -1, threads = 1, sampleChunks = -1, sample = -1;
        bool ss = false;
        enum { REF_SEQUENCE = 1000, START, LENGTH, ROOT_GENOME, EXCLUDE_GENOMES, LEAF_INSERTIONS, NO_ANCESTOR_BASES, CHUNK, THREADS, SS,
               SAMPLE_CHUNKS, SAMPLE };
        static const option longOpts[] = {{"refSequence", required_argument, nullptr, REF_SEQUENCE},
                                          {"start", required_argument, nullptr, START},
                                          {"length", required_argument, nullptr, LENGTH},
                                          {"rootGenome", required_argument, nullptr, ROOT_GENOME},
                                          {"excludeGenomes", required_argument, nullptr, EXCLUDE_GENOMES},
                                          {"leafInsertions", no_argument, nullptr, LEAF_INSERTIONS},
                                          {"noAncestorBases", no_argument, nullptr, NO_ANCESTOR_BASES},
                                          {"chunk", required_argument, nullptr, CHUNK},
                                          {"threads", required_argument, nullptr, THREADS},
                                          {"ss", no_argument, nullptr, SS},
                                          {"sampleChunks", required_argument, nullptr, SAMPLE_CHUNKS},
                                          {"sample", required_argument, nullptr, SAMPLE},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        optind = 1;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case REF_SEQUENCE: refSeqName = optarg; break;
            case START: start = parseInt("start", optarg); break;
            case LENGTH: length = parseInt("length", optarg); break;
            case ROOT_GENOME: opts.rootGenome = optarg; break;
            case EXCLUDE_GENOMES: {
                stringstream ss(optarg);
                string x;
                while (getline(ss, x, ','))
                    if (!x.empty()) opts.exclude.push_back(x);
                break;
            }
            case LEAF_INSERTIONS: opts.leafInsertions = true; break;
            case NO_ANCESTOR_BASES: opts.ancestorBases = false; break;
            case CHUNK: opts.chunk = parseInt("chunk", optarg); break;
            case THREADS: threads = parseInt("threads", optarg); break;
            case SS: ss = true; break;
            case SAMPLE_CHUNKS: sampleChunks = parseInt("sampleChunks", optarg); break;
            case SAMPLE: sample = parseInt("sample", optarg); break;
            case 'h': cout << COLUMNS_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        if (threads < 1 || opts.chunk < 1) throw CliError("--threads and --chunk must be positive");
        if (!ss && (sampleChunks != -1 || sample != -1)) throw CliError("--sampleChunks and --sample need --ss");
        if ((sampleChunks != -1 && sampleChunks < 1) || (sample != -1 && sample < 1)) throw CliError("--sampleChunks and --sample must be positive");
        if (ss) opts.ancestorBases = false;   /* only leaves are rows */
        Alignment aln(argv[optind]);
        const string outPath = argv[optind + 1];
        const Genome *backbone = opts.rootGenome.empty() ? aln.root() : aln.genome(opts.rootGenome);
        if (!backbone) throw runtime_error("no genome " + opts.rootGenome);
        vector<RefInterval> ranges;
        if (!refSeqName.empty()) {
            int64_t i = backbone->findSequence(refSeqName);
            if (i < 0) throw runtime_error("no sequence " + refSeqName + " in " + backbone->name());
            const int64_t len = backbone->sequenceLength((uint32_t)i);
            if (start < 0 || start > len) throw CliError("--start is outside " + refSeqName);
            const int64_t e = length < 0 ? len : min(len, start + length);
            if (start < e) ranges.push_back(RefInterval{(uint32_t)i, start, e});
            if (ranges.empty()) throw CliError("the range is empty");
        } else if (start != 0 || length >= 0) {
            throw CliError("--start / --length need --refSequence");
        }
        UniversalColumns cols(aln, ranges, opts, (int)threads);
        if (ss) {
            /* the leaves' rows, in preorder (as hal 4d's), by phast name */
            vector<int32_t> rowOf(cols.genomes(), -1);
            vector<string> rows;
            map<string, string> seenNames;
            for (uint32_t i = 0; i < cols.genomes(); i++) {
                const Genome *g = cols.genome(i);
                if (!g->children().empty()) continue;
                string pn = g->name().substr(0, g->name().find('.'));
                if (seenNames.count(pn)) throw runtime_error("leaves " + seenNames[pn] + " and " + g->name() + " have the same phast name " + pn);
                seenNames[pn] = g->name();
                rowOf[i] = (int32_t)rows.size();
                rows.push_back(pn);
            }
            vector<size_t> which;
            const size_t nc = cols.chunks();
            if (sampleChunks < 0 || (size_t)sampleChunks >= nc) {
                for (size_t k = 0; k < nc; k++) which.push_back(k);
            } else {
                for (int64_t j = 0; j < sampleChunks; j++) which.push_back((size_t)((double)j * (double)nc / (double)sampleChunks));
                which.erase(unique(which.begin(), which.end()), which.end());
            }
            struct Site {
                uint64_t hash;
                uint32_t genome, seq;
                int64_t pos;
                string tuple;
            };
            vector<Site> sites;
            uint64_t parts = 0, alone = 0;
            vector<unique_ptr<ColumnChunk>> cks((size_t)threads);
            orderedParallel<vector<Site>>(
                which.size(), (int)threads,
                [&](size_t j, int w) {
                    unique_ptr<ColumnChunk> &ch = cks[(size_t)w];
                    if (!ch) ch.reset(new ColumnChunk());
                    cols.chunk(which[j], w, *ch);
                    vector<Site> out;
                    vector<int32_t> root;
                    map<pair<int32_t, uint32_t>, int> kids;
                    map<int32_t, vector<int64_t>> part;
                    for (int64_t u = 0; u < ch->cols; u++) {
                        const int64_t e0 = ch->colStart[(size_t)u], e1 = ch->colStart[(size_t)u + 1];
                        kids.clear();
                        part.clear();
                        root.assign((size_t)(e1 - e0), -1);
                        for (int64_t e = e0; e < e1; e++)
                            if (ch->parent[(size_t)e] >= 0) kids[make_pair(ch->parent[(size_t)e], ch->genome[(size_t)e])]++;
                        for (int64_t e = e0; e < e1; e++) {
                            const int32_t i = (int32_t)(e - e0), q = ch->parent[(size_t)e];
                            root[(size_t)i] = q < 0 || kids[make_pair(q, ch->genome[(size_t)e])] > 1 ? i : root[(size_t)q];
                            if (rowOf[ch->genome[(size_t)e]] >= 0 && tupleBase(ch->base[(size_t)e]) != 'N') part[root[(size_t)i]].push_back(e);
                        }
                        for (auto &pt : part) {
                            if (pt.second.size() < 2) continue;
                            const int64_t r = e0 + pt.first;
                            const ColumnPiece &p = ch->pieces[ch->piece[(size_t)r]];
                            const int64_t d = u - p.u, pos = p.reversed ? p.start + p.length - 1 - d : p.start + d;
                            const Genome *g = cols.genome(p.genome);
                            Site s{fnvHash(g->name() + '\0' + g->sequenceName(p.seq) + '\0' + to_string(pos)), p.genome, p.seq, pos,
                                   string(rows.size(), '*')};
                            for (int64_t e : pt.second) s.tuple[(size_t)rowOf[ch->genome[(size_t)e]]] = tupleBase(ch->base[(size_t)e]);
                            out.push_back(move(s));
                        }
                    }
                    return out;
                },
                [&](size_t, vector<Site> &v) {
                    parts += v.size();
                    for (Site &x : v) sites.push_back(move(x));
                });
            (void)alone;
            sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) { return tie(a.hash, a.genome, a.seq, a.pos) < tie(b.hash, b.genome, b.seq, b.pos); });
            if (sample > 0 && (size_t)sample < sites.size()) sites.resize((size_t)sample);
            sort(sites.begin(), sites.end(), [](const Site &a, const Site &b) { return tie(a.genome, a.seq, a.pos) < tie(b.genome, b.seq, b.pos); });
            vector<pair<const string *, uint64_t>> tuples;
            unordered_map<string, size_t> tupleOf;
            for (const Site &x : sites) {
                auto it = tupleOf.find(x.tuple);
                if (it == tupleOf.end()) {
                    tupleOf[x.tuple] = tuples.size();
                    tuples.push_back(make_pair(&x.tuple, 1));
                } else {
                    tuples[it->second].second++;
                }
            }
            ofstream o(outPath);
            if (!o) throw runtime_error("cannot write " + outPath);
            o << "NSEQS = " << rows.size() << "\nLENGTH = " << sites.size() << "\nTUPLE_SIZE = 1\nNTUPLES = " << tuples.size() << "\nNAMES = ";
            for (size_t i = 0; i < rows.size(); i++) o << (i ? "," : "") << rows[i];
            o << "\nALPHABET = ACGT\nNCATS = -1\n\n";
            for (size_t i = 0; i < tuples.size(); i++) o << i << '\t' << *tuples[i].first << '\t' << tuples[i].second << '\n';
            o.close();
            if (!o) throw runtime_error("cannot write " + outPath);
            cerr << "hal columns: " << which.size() << " of " << nc << " chunks, " << parts << " single-copy parts with two leaves or more, "
                 << sites.size() << " written (" << tuples.size() << " distinct)" << endl;
            return 0;
        }
        vector<string> names(cols.genomes());
        for (uint32_t i = 0; i < cols.genomes(); i++) names[i] = cols.genome(i)->name() + ".";
        ofstream file;
        if (outPath != "stdout") {
            file.open(outPath);
            if (!file) throw runtime_error("cannot write " + outPath);
        }
        ostream &out = outPath == "stdout" ? cout : file;
        vector<unique_ptr<ColumnChunk>> chunks((size_t)threads);
        orderedParallel<string>(
            cols.chunks(), (int)threads,
            [&](size_t k, int w) {
                unique_ptr<ColumnChunk> &ch = chunks[(size_t)w];
                if (!ch) ch.reset(new ColumnChunk());
                cols.chunk(k, w, *ch);
                ostringstream s;
                for (int64_t u = 0; u < ch->cols; u++) {
                    s << *ch->axisName << '\t' << ch->axisStart + u;
                    const int64_t e0 = ch->colStart[(size_t)u], e1 = ch->colStart[(size_t)u + 1];
                    for (int64_t e = e0; e < e1; e++) {
                        const ColumnPiece &p = ch->pieces[ch->piece[(size_t)e]];
                        const Genome *g = cols.genome(p.genome);
                        const int64_t d = u - p.u, pos = p.reversed ? p.start + p.length - 1 - d : p.start + d;
                        s << '\t' << names[p.genome] << g->sequenceName(p.seq) << '|' << pos << '|' << (p.reversed ? '-' : '+') << '|'
                          << ch->base[(size_t)e] << '|';
                        if (ch->parent[(size_t)e] < 0) s << '-';
                        else s << ch->parent[(size_t)e];
                        s << '|' << ch->piece[(size_t)e];
                    }
                    s << '\n';
                }
                return s.str();
            },
            [&](size_t, string &text) { out << text; });
        out.flush();
        if (!out) throw runtime_error("cannot write " + outPath);
        cerr << "hal columns: " << cols.chunks() << " chunks, " << cols.columns() << " columns" << endl;
        return 0;
    });
}
}
