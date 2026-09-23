/*
 * hal3ToMaf: MAF from a HAL format 3 alignment, row by row (see hal3MafExport.h).
 */
#include "hal3MafExport.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

static void usage() {
    cerr << "usage: hal3ToMaf alignment_dir out.maf [options]\n"
            "MAF blocks anchored on a reference genome, from a HAL format 3 alignment.  out.maf may be\n"
            "stdout.  Written BGZF-compressed to a file, plain to stdout, unless told otherwise.\n"
            "  --refGenome G        reference genome (default: the root)\n"
            "  --refSequence S      only this reference sequence\n"
            "  --start N --length N only this range of it (0-based)\n"
            "  --refTargets BED     only these reference intervals (BED3+)\n"
            "  --noAncestors        no rows for ancestral genomes\n"
            "  --noDupes            one copy per base (the canonical paralog)\n"
            "  --onlyOrthologs      rows from liftover to each genome's common ancestor with the reference\n"
            "  --targetGenomes a,b  only rows for these genomes (and the reference)\n"
            "  --rootGenome R       only rows for the subtree of R (and the reference)\n"
            "  --breakFrac F        cut a block where this share of rows break (default 0.5; 0 = at every break)\n"
            "  --gapBudget G        cut a block before gap cells exceed G x aligned bases (default 1.0)\n"
            "  --maxBlockCols N     at most N columns per block (default 10000)\n"
            "  --maxBlockCells N    at most N cells (rows x columns) per block (default 5000000)\n"
            "  --noDna              write N for every base\n"
            "  --onlySequenceNames  name rows by sequence alone, not genome.sequence\n"
            "  --uncompressed       plain text to a file\n"
            "  --bgzip              BGZF to stdout\n";
    exit(1);
}

int main(int argc, char **argv) {
    vector<string> args;
    hal3::MafOptions opts;
    string refName, refSeqName, bedPath;
    int64_t start = 0, length = -1;
    bool uncompressed = false, bgzipStdout = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        auto val = [&]() -> string {
            if (i + 1 >= argc) usage();
            return argv[++i];
        };
        if (a == "--refGenome") refName = val();
        else if (a == "--refSequence") refSeqName = val();
        else if (a == "--start") start = atoll(val().c_str());
        else if (a == "--length") length = atoll(val().c_str());
        else if (a == "--refTargets") bedPath = val();
        else if (a == "--noAncestors") opts.ancestors = false;
        else if (a == "--noDupes") opts.dupes = false;
        else if (a == "--onlyOrthologs") opts.onlyOrthologs = true;
        else if (a == "--targetGenomes") {
            stringstream ss(val());
            string t;
            while (getline(ss, t, ',')) opts.targets.push_back(t);
        } else if (a == "--rootGenome") opts.rootGenome = val();
        else if (a == "--breakFrac") opts.breakFrac = atof(val().c_str());
        else if (a == "--gapBudget") opts.gapBudget = atof(val().c_str());
        else if (a == "--maxBlockCols") opts.maxBlockCols = atoll(val().c_str());
        else if (a == "--maxBlockCells") opts.maxBlockCells = atoll(val().c_str());
        else if (a == "--chunk") opts.chunk = atoll(val().c_str());
        else if (a == "--noDna") opts.noDna = true;
        else if (a == "--onlySequenceNames") opts.onlySequenceNames = true;
        else if (a == "--uncompressed") uncompressed = true;
        else if (a == "--bgzip") bgzipStdout = true;
        else if (a == "--maxRefGap" || a == "--maxBlockLen") {
            cerr << "hal3ToMaf: " << a << " is not supported; see --maxBlockCols / --maxBlockCells ("
                 << "insertion columns are not implemented yet)" << endl;
            return 1;
        } else if (a.size() > 1 && a[0] == '-' && a[1] == '-') usage();
        else args.push_back(a);
    }
    if (args.size() != 2 || opts.maxBlockCols < 1 || opts.maxBlockCells < 1 || opts.chunk < 1) usage();
    try {
        hal3::Alignment aln(args[0]);
        const hal3::Genome *ref = refName.empty() ? aln.root() : aln.genome(refName);
        if (!ref) throw runtime_error("no genome " + refName);
        bool toStdout = args[1] == "stdout";
        unique_ptr<hal3::MafSink> sink = hal3::openMafSink(args[1], toStdout ? bgzipStdout : !uncompressed);
        hal3::MafExporter exporter(aln, ref, opts, *sink);
        exporter.writeHeader();
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
                if (line.empty() || line[0] == '#' || !line.compare(0, 5, "track") || !line.compare(0, 7, "browser")) continue;
                stringstream ss(line);
                string chrom;
                int64_t s, e;
                if (!(ss >> chrom >> s >> e)) throw runtime_error("bad BED line: " + line);
                uint32_t i = seqIndex(chrom);
                e = min<int64_t>(e, ref->sequenceLength(i));
                if (s < e) exporter.exportInterval(i, s, e);
            }
        } else if (!refSeqName.empty()) {
            uint32_t i = seqIndex(refSeqName);
            int64_t len = ref->sequenceLength(i);
            int64_t e = length < 0 ? len : min(len, start + length);
            if (start < 0 || start > len) throw runtime_error("--start is outside " + refSeqName);
            exporter.exportInterval(i, start, e);
        } else {
            if (start != 0 || length >= 0) throw runtime_error("--start / --length need --refSequence");
            for (uint32_t i = 0; i < ref->numSequences(); i++) exporter.exportInterval(i, 0, ref->sequenceLength(i));
        }
        sink->close();
        const hal3::MafExporter::Stats &st = exporter.stats();
        cerr << "hal3ToMaf: " << st.blocks << " blocks, " << st.columns << " columns, " << st.rows << " rows, "
             << st.alignedBases << " aligned bases" << endl;
    } catch (exception &e) {
        cerr << "hal3ToMaf: " << e.what() << endl;
        return 1;
    }
    return 0;
}
