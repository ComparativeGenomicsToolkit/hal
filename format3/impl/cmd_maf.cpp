/*
 * hal maf: MAF from an alignment, row by row (see hal3MafExport.h).  Positional
 * arguments as hal2maf.
 */
#include "hal3Cli.h"
#include "hal3MafExport.h"

#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

namespace hal3 {

static const char *MAF_USAGE =
    "usage: hal maf [options] <alignment> <out.maf>\n"
    "MAF (or TAF) blocks anchored on a reference genome.  out.maf may be stdout.  The name picks\n"
    "the format: .taf or .taf.gz is TAF, anything else MAF; .gz is BGZF-compressed (readable by\n"
    "gzip, indexable by taffy index), anything else plain.  Every column holding a reference\n"
    "base is hal2maf's; blocks are cut differently (see --breakFrac, --gapBudget) and joined\n"
    "across insertions (--joinGap), as taffy norm does, with homologous inserted bases aligned.\n"
    "Rows: the reference, then genomes by tree proximity to it (as cactus-hal2maf orders them).\n"
    "The output is the same for any number of threads.\n"
    "options:\n"
    "  --refGenome G        reference genome (default: the root)\n"
    "  --refSequence S      only this reference sequence\n"
    "  --start N            only from this position of it (0-based; with --refSequence)\n"
    "  --length N           only this many bases (with --refSequence)\n"
    "  --refTargets BED     only these reference intervals (BED3+), in file order (for output\n"
    "                       taffy index can use, give each sequence's intervals together)\n"
    "  --noAncestors        no rows for ancestral genomes\n"
    "  --noDupes            one copy per base (the canonical paralog)\n"
    "  --onlyOrthologs      rows from liftover to each genome's common ancestor with the reference\n"
    "  --targetGenomes a,b  only rows for these genomes (and the reference)\n"
    "  --rootGenome R       only rows for the subtree of R (and the reference)\n"
    "  --excludeGenomes a,b no rows for these genomes (e.g. _MINIGRAPH_).  The walk still goes\n"
    "                       through them: other rows' reference columns are unchanged, but\n"
    "                       blocks and insertion columns are formed without them\n"
    "  --refSingleCopy      no rows for the reference genome's other copies of its bases\n"
    "  --joinGap N          continue a row across up to N inserted bases (bases in no column\n"
    "                       with a reference base), which get columns of their own, shared by\n"
    "                       homologous bases (default 30; 0: no insertion columns)\n"
    "  --breakFrac F        cut a block where this share of its rows break (default 0.5; 0 cuts at\n"
    "                       every break)\n"
    "  --gapBudget G        cut a block before its gap cells exceed G x its aligned bases (default 1)\n"
    "  --maxBlockCols N     at most N columns per block (default 10000)\n"
    "  --maxBlockCells N    at most N cells, rows x columns, per block (default 5000000)\n"
    "  --noDna              write N for every base\n"
    "  --onlySequenceNames  name rows by sequence alone, not genome.sequence\n"
    "  --threads N          export N chunks of the reference at once (default 1)\n"
    "  --taf                write TAF, whatever the name (e.g. to stdout)\n"
    "  --noRunLength        TAF: write each column's bases in full, not run-length encoded\n"
    "  --bgzip              BGZF, whatever the name (e.g. to stdout)\n"
    "  --uncompressed       plain text, whatever the name\n"
    "  -h, --help\n";

int main_maf(int argc, char **argv) {
    return runCommand("maf", [&]() {
        MafOptions opts;
        MafFormat format;
        string refName, refSeqName, bedPath;
        int64_t start = 0, length = -1, threads = 1;
        bool uncompressed = false, bgzip = false, taf = false;
        enum {
            REF_GENOME = 1000, REF_SEQUENCE, START, LENGTH, REF_TARGETS, NO_ANCESTORS, NO_DUPES, ONLY_ORTHOLOGS,
            TARGET_GENOMES, ROOT_GENOME, BREAK_FRAC, GAP_BUDGET, MAX_BLOCK_COLS, MAX_BLOCK_CELLS, CHUNK, NO_DNA,
            ONLY_SEQUENCE_NAMES, UNCOMPRESSED, BGZIP, MAX_REF_GAP, MAX_BLOCK_LEN, JOIN_GAP, THREADS, TAF,
            NO_RUN_LENGTH, EXCLUDE_GENOMES, REF_SINGLE_COPY
        };
        static const option longOpts[] = {
            {"refGenome", required_argument, nullptr, REF_GENOME},
            {"refSequence", required_argument, nullptr, REF_SEQUENCE},
            {"start", required_argument, nullptr, START},
            {"length", required_argument, nullptr, LENGTH},
            {"refTargets", required_argument, nullptr, REF_TARGETS},
            {"noAncestors", no_argument, nullptr, NO_ANCESTORS},
            {"noDupes", no_argument, nullptr, NO_DUPES},
            {"onlyOrthologs", no_argument, nullptr, ONLY_ORTHOLOGS},
            {"targetGenomes", required_argument, nullptr, TARGET_GENOMES},
            {"rootGenome", required_argument, nullptr, ROOT_GENOME},
            {"breakFrac", required_argument, nullptr, BREAK_FRAC},
            {"gapBudget", required_argument, nullptr, GAP_BUDGET},
            {"maxBlockCols", required_argument, nullptr, MAX_BLOCK_COLS},
            {"maxBlockCells", required_argument, nullptr, MAX_BLOCK_CELLS},
            {"chunk", required_argument, nullptr, CHUNK},
            {"noDna", no_argument, nullptr, NO_DNA},
            {"onlySequenceNames", no_argument, nullptr, ONLY_SEQUENCE_NAMES},
            {"uncompressed", no_argument, nullptr, UNCOMPRESSED},
            {"bgzip", no_argument, nullptr, BGZIP},
            {"joinGap", required_argument, nullptr, JOIN_GAP},
            {"threads", required_argument, nullptr, THREADS},
            {"taf", no_argument, nullptr, TAF},
            {"noRunLength", no_argument, nullptr, NO_RUN_LENGTH},
            {"excludeGenomes", required_argument, nullptr, EXCLUDE_GENOMES},
            {"refSingleCopy", no_argument, nullptr, REF_SINGLE_COPY},
            {"maxRefGap", required_argument, nullptr, MAX_REF_GAP},
            {"maxBlockLen", required_argument, nullptr, MAX_BLOCK_LEN},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case REF_GENOME: refName = optarg; break;
            case REF_SEQUENCE: refSeqName = optarg; break;
            case START: start = parseInt("start", optarg); break;
            case LENGTH: length = parseInt("length", optarg); break;
            case REF_TARGETS: bedPath = optarg; break;
            case NO_ANCESTORS: opts.ancestors = false; break;
            case NO_DUPES: opts.dupes = false; break;
            case ONLY_ORTHOLOGS: opts.onlyOrthologs = true; break;
            case TARGET_GENOMES:
            case EXCLUDE_GENOMES: {
                stringstream ss(optarg);
                string t;
                while (getline(ss, t, ','))
                    if (!t.empty()) (c == TARGET_GENOMES ? opts.targets : opts.exclude).push_back(t);
                break;
            }
            case REF_SINGLE_COPY: opts.refSingleCopy = true; break;
            case ROOT_GENOME: opts.rootGenome = optarg; break;
            case BREAK_FRAC: opts.breakFrac = parseDouble("breakFrac", optarg); break;
            case GAP_BUDGET: opts.gapBudget = parseDouble("gapBudget", optarg); break;
            case MAX_BLOCK_COLS: opts.maxBlockCols = parseInt("maxBlockCols", optarg); break;
            case MAX_BLOCK_CELLS: opts.maxBlockCells = parseInt("maxBlockCells", optarg); break;
            case CHUNK: opts.chunk = parseInt("chunk", optarg); break;
            case NO_DNA: opts.noDna = true; break;
            case ONLY_SEQUENCE_NAMES: opts.onlySequenceNames = true; break;
            case UNCOMPRESSED: uncompressed = true; break;
            case BGZIP: bgzip = true; break;
            case JOIN_GAP: opts.joinGap = parseInt("joinGap", optarg); break;
            case THREADS: threads = parseInt("threads", optarg); break;
            case TAF: taf = true; break;
            case NO_RUN_LENGTH: format.runLength = false; break;
            case MAX_REF_GAP:
            case MAX_BLOCK_LEN:
                throw CliError(c == MAX_REF_GAP ? "--maxRefGap is not supported: see --joinGap"
                                                : "--maxBlockLen is not supported: see --maxBlockCols / --maxBlockCells");
            case 'h': cout << MAF_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        if (opts.maxBlockCols < 1 || opts.maxBlockCells < 1 || opts.chunk < 1 || opts.joinGap < 0 || threads < 1)
            throw CliError("--maxBlockCols, --maxBlockCells, --chunk and --threads must be positive, --joinGap not negative");
        if (bgzip && uncompressed) throw CliError("--bgzip and --uncompressed cannot be combined");
        string dir = argv[optind], outPath = argv[optind + 1];
        auto endsWith = [](const string &s, const string &t) { return s.size() >= t.size() && !s.compare(s.size() - t.size(), t.size(), t); };
        bool gz = endsWith(outPath, ".gz");
        format.bgzf = bgzip || (gz && !uncompressed);
        format.taf = taf || endsWith(gz ? outPath.substr(0, outPath.size() - 3) : outPath, ".taf");
        Alignment aln(dir);
        const Genome *ref = refName.empty() ? aln.root() : aln.genome(refName);
        if (!ref) throw runtime_error("no genome " + refName);
        vector<RefInterval> intervals;
        auto seqIndex = [&](const string &name) {
            int64_t i = ref->findSequence(name);
            if (i < 0) throw runtime_error("no sequence " + name + " in " + ref->name());
            return (uint32_t)i;
        };
        if (!bedPath.empty()) {
            ifstream bed(bedPath);
            if (!bed) throw runtime_error("cannot open " + bedPath);
            string line;
            while (getline(bed, line)) {
                if (line.empty() || line[0] == '#' || !line.compare(0, 5, "track") || !line.compare(0, 7, "browser"))
                    continue;
                stringstream ss(line);
                string chrom;
                int64_t s, e;
                if (!(ss >> chrom >> s >> e)) throw runtime_error("bad BED line: " + line);
                uint32_t i = seqIndex(chrom);
                e = min<int64_t>(e, ref->sequenceLength(i));
                if (s < e) intervals.push_back(RefInterval{i, s, e});
            }
        } else if (!refSeqName.empty()) {
            uint32_t i = seqIndex(refSeqName);
            int64_t len = ref->sequenceLength(i);
            int64_t e = length < 0 ? len : min(len, start + length);
            if (start < 0 || start > len) throw runtime_error("--start is outside " + refSeqName);
            intervals.push_back(RefInterval{i, start, e});
        } else {
            if (start != 0 || length >= 0) throw CliError("--start / --length need --refSequence");
            for (uint32_t i = 0; i < ref->numSequences(); i++) intervals.push_back(RefInterval{i, 0, ref->sequenceLength(i)});
        }
        MafStats st = exportMaf(aln, ref, opts, format, intervals, outPath, (int)threads);
        cerr << "hal maf: " << st.blocks << " blocks, " << st.columns << " columns, " << st.rows << " rows, "
             << st.alignedBases << " aligned bases; " << st.insertionColumns << " insertion columns ("
             << st.insertedBases << " bases), " << st.joinsRejected << " joins refused" << endl;
        return 0;
    });
}
}
