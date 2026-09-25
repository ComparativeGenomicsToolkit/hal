/*
 * hal fasta: genome sequences as FASTA, as hal2fasta (same positional arguments and
 * output), read straight from the genomes' .2bit files.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "twobit64.h"

#include <cctype>
#include <cstdio>
#include <cstring>
#include <deque>
#include <getopt.h>
#include <iostream>
#include <memory>

using namespace std;

namespace hal3 {

static const char *FASTA_USAGE =
    "usage: hal fasta [options] <alignment> <genome>\n"
    "The sequences of a genome (or a subtree of genomes) as FASTA, as hal2fasta.\n"
    "options:\n"
    "  --outFaPath F         output file (default stdout)\n"
    "  --sequence S          only this sequence\n"
    "  --start N             from this position (of the sequence, or of the genome's\n"
    "                        concatenated sequences)\n"
    "  --length N            this many bases (default: to the end)\n"
    "  --lineWidth N         bases per line (default 80)\n"
    "  --subtree             every genome in the subtree of <genome>, breadth first\n"
    "  --upper               upper-case every base\n"
    "  --ucscSequenceNames   name sequences genome.sequence\n"
    "  -h, --help\n";

namespace {

struct Writer {
    FILE *out;
    int64_t lineWidth;
    bool upper, fullNames;
    string buf;

    /* One sequence's [start, start+length) (length 0: to the end), as hal2fasta's
     * printSequence. */
    void sequence(const Genome *g, tb_file *dna, int64_t tbIndex, uint32_t seq, int64_t start, int64_t length) {
        int64_t seqLen = g->sequenceLength(seq);
        if (start > seqLen)
            throw runtime_error("Start position " + to_string(start) + " is out of range for sequence " +
                                g->sequenceName(seq) + ", which has length " + to_string(seqLen));
        if (length == 0) length = seqLen - start;
        int64_t last = start + length;
        if (last > seqLen)
            throw runtime_error("Specified range [" + to_string(start) + "," + to_string(length) +
                                "] is out of range for sequence " + g->sequenceName(seq) + ", which has length " +
                                to_string(seqLen));
        fprintf(out, ">%s%s%s\n", fullNames ? g->name().c_str() : "", fullNames ? "." : "", g->sequenceName(seq));
        const int64_t piece = (int64_t)1 << 20;
        for (int64_t p = start; p < last; p += piece) {
            int64_t e = min(last, p + piece);
            buf.resize((size_t)(e - p));
            if (tb_fetch(dna, (uint32_t)tbIndex, p, e, &buf[0], 0) != 0)
                throw runtime_error(string("cannot read DNA of ") + g->name() + "." + g->sequenceName(seq));
            if (upper)
                for (char &ch : buf) ch = (char)toupper((unsigned char)ch);
            /* lines of lineWidth counted from start */
            for (int64_t q = p; q < e;) {
                int64_t lineEnd = min(e, start + ((q - start) / lineWidth + 1) * lineWidth);
                fwrite(&buf[q - p], 1, (size_t)(lineEnd - q), out);
                if (lineEnd == start + ((q - start) / lineWidth + 1) * lineWidth || lineEnd == last) fputc('\n', out);
                q = lineEnd;
            }
        }
    }

    /* A genome, or its range of concatenated sequences, as hal2fasta's printGenome. */
    void genome(const Alignment &aln, const Genome *g, int64_t seqOnly, int64_t start, int64_t length) {
        const GenomeEntry *e = aln.manifest().find(g->name());
        unique_ptr<tb_file, void (*)(tb_file *)> dna(aln.openDna(g), tb_close);
        auto tbIndex = [&](uint32_t i) {
            int64_t t = i < tb_n_seqs(dna.get()) && !strcmp(tb_name(dna.get(), i), g->sequenceName(i))
                            ? (int64_t)i
                            : tb_find(dna.get(), g->sequenceName(i));
            if (t < 0) throw runtime_error(string("sequence ") + g->sequenceName(i) + " is not in " + e->dnaPath);
            return t;
        };
        if (seqOnly >= 0) {
            sequence(g, dna.get(), tbIndex((uint32_t)seqOnly), (uint32_t)seqOnly, start, length);
            return;
        }
        int64_t total = g->sequenceStart(g->numSequences() - 1) + g->sequenceLength(g->numSequences() - 1);
        if (start + length > total)
            throw runtime_error("Specified range [" + to_string(start) + "," + to_string(length) +
                                "] is out of range for genome " + g->name() + ", which has length " + to_string(total));
        if (length == 0) length = total - start;
        for (uint32_t i = 0; i < g->numSequences(); i++) {
            int64_t s0 = g->sequenceStart(i), len = g->sequenceLength(i);
            if (start + length > s0 && start < s0 + len) {
                int64_t readStart = start < s0 ? 0 : start - s0;
                int64_t readLen = start + length >= s0 + len ? len - readStart : start + length - s0 - readStart;
                sequence(g, dna.get(), tbIndex(i), i, readStart, readLen);
            } else if (start + length <= s0) {
                break;
            }
        }
    }
};
}

int main_fasta(int argc, char **argv) {
    return runCommand("fasta", [&]() {
        string outPath = "stdout", seqName;
        int64_t start = 0, length = 0, lineWidth = 80;
        bool subtree = false, upper = false, fullNames = false, startSet = false, lengthSet = false;
        enum { OUT = 1000, SEQUENCE, START, LENGTH, LINE_WIDTH, SUBTREE, UPPER, UCSC };
        static const option longOpts[] = {{"outFaPath", required_argument, nullptr, OUT},
                                          {"sequence", required_argument, nullptr, SEQUENCE},
                                          {"start", required_argument, nullptr, START},
                                          {"length", required_argument, nullptr, LENGTH},
                                          {"lineWidth", required_argument, nullptr, LINE_WIDTH},
                                          {"subtree", no_argument, nullptr, SUBTREE},
                                          {"upper", no_argument, nullptr, UPPER},
                                          {"ucscSequenceNames", no_argument, nullptr, UCSC},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case OUT: outPath = optarg; break;
            case SEQUENCE: seqName = optarg; break;
            case START: start = parseInt("start", optarg); startSet = true; break;
            case LENGTH: length = parseInt("length", optarg); lengthSet = true; break;
            case LINE_WIDTH: lineWidth = parseInt("lineWidth", optarg); break;
            case SUBTREE: subtree = true; break;
            case UPPER: upper = true; break;
            case UCSC: fullNames = true; break;
            case 'h': cout << FASTA_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        if (lineWidth < 1 || start < 0 || length < 0) throw CliError("--lineWidth must be positive, --start and --length not negative");
        if (subtree && (startSet || lengthSet || !seqName.empty()))
            throw CliError("--start, --length and --sequence cannot be used with --subtree");
        Alignment aln(argv[optind]);
        const Genome *g = aln.genome(argv[optind + 1]);
        if (!g) throw runtime_error(string("Genome ") + argv[optind + 1] + " not found");
        int64_t seqOnly = -1;
        if (!seqName.empty()) {
            seqOnly = g->findSequence(seqName);
            if (seqOnly < 0) throw runtime_error("Sequence " + seqName + " not found");
        }
        FILE *out = outPath == "stdout" ? stdout : fopen(outPath.c_str(), "w");
        if (!out) throw runtime_error("cannot write " + outPath);
        Writer w{out, lineWidth, upper, fullNames, {}};
        deque<const Genome *> queue{g};
        while (!queue.empty()) {
            const Genome *cur = queue.front();
            queue.pop_front();
            w.genome(aln, cur, seqOnly, start, length);
            if (subtree)
                for (const Genome *ch : cur->children()) queue.push_back(ch);
        }
        if ((out == stdout ? fflush(out) : fclose(out)) != 0) throw runtime_error("write to " + outPath + " failed");
        return 0;
    });
}
}
