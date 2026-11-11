/*
 * hal columns: universal columns as trees of bases (hal3Columns.h), as text -- to check
 * them, and to look at them.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Columns.h"

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <sstream>

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
    "  --threads N          chunks on N threads (default 1); the output is the same\n";

int main_columns(int argc, char **argv) {
    return runCommand("columns", [&]() {
        ColumnOptions opts;
        string refSeqName;
        int64_t start = 0, length = -1, threads = 1;
        enum { REF_SEQUENCE = 1000, START, LENGTH, ROOT_GENOME, EXCLUDE_GENOMES, LEAF_INSERTIONS, NO_ANCESTOR_BASES, CHUNK, THREADS };
        static const option longOpts[] = {{"refSequence", required_argument, nullptr, REF_SEQUENCE},
                                          {"start", required_argument, nullptr, START},
                                          {"length", required_argument, nullptr, LENGTH},
                                          {"rootGenome", required_argument, nullptr, ROOT_GENOME},
                                          {"excludeGenomes", required_argument, nullptr, EXCLUDE_GENOMES},
                                          {"leafInsertions", no_argument, nullptr, LEAF_INSERTIONS},
                                          {"noAncestorBases", no_argument, nullptr, NO_ANCESTOR_BASES},
                                          {"chunk", required_argument, nullptr, CHUNK},
                                          {"threads", required_argument, nullptr, THREADS},
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
            case 'h': cout << COLUMNS_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        if (threads < 1 || opts.chunk < 1) throw CliError("--threads and --chunk must be positive");
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
