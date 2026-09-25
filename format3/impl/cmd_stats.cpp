/*
 * hal stats: facts about an alignment, as halStats (same options and output where the
 * idea carries over).  Format 3 has no segments: where halStats counts top and bottom
 * segments, this counts runs -- maximal gapless aligned blocks of a genome's branch to
 * its parent.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "twobit64.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>

using namespace std;

namespace hal3 {

static const char *STATS_USAGE =
    "usage: hal stats [option] <alignment>\n"
    "Facts about an alignment, as halStats.  With no option, a summary: the tree, and per\n"
    "genome its children, length, sequences, runs (gapless aligned blocks of its branch to\n"
    "its parent) and file sizes.  At most one option:\n"
    "  --genomes              the genomes\n"
    "  --tree                 the Newick tree\n"
    "  --root                 the root genome\n"
    "  --children G           G's children\n"
    "  --parent G             G's parent\n"
    "  --branchLength G       the length of the branch above G\n"
    "  --branches             the branches, named by their child genomes\n"
    "  --span A,B,...         the branches of the spanning tree of these genomes: the path\n"
    "                         from A to B for two, tree order for more\n"
    "  --spanRoot A,B,...     as --span, with the spanning tree's root\n"
    "  --sequences G          G's sequences, comma separated\n"
    "  --sequenceStats G      per sequence of G: length, runs, bases aligned to the parent\n"
    "  --chromSizes G         per sequence of G: name and length\n"
    "  --bedSequences G       G's sequences as BED\n"
    "  --runs G               the runs of the branch above G, as BED: sequence, start, end,\n"
    "                         parentSequence:parentStart, 0, strand\n"
    "  --baseComp G,step      base composition of G, sampling every step bases\n"
    "  --genomeMetaData G     G's metadata\n"
    "  --metaData             the alignment's metadata\n"
    "  --assembly             assembly statistics for every genome\n"
    "  -h, --help\n";

namespace {

vector<const Genome *> preorder(const Alignment &aln) {
    vector<const Genome *> out, stack{aln.root()};
    while (!stack.empty()) {
        const Genome *g = stack.back();
        stack.pop_back();
        out.push_back(g);
        for (auto c = g->children().rbegin(); c != g->children().rend(); ++c) stack.push_back(*c);
    }
    return out;
}

const Genome *need(const Alignment &aln, const string &name) {
    const Genome *g = aln.genome(name);
    if (!g) throw runtime_error("Genome " + name + " not found.");
    return g;
}

vector<string> split(const string &s, char sep) {
    vector<string> out;
    stringstream ss(s);
    string t;
    while (getline(ss, t, sep))
        if (!t.empty()) out.push_back(t);
    return out;
}

/* Runs and aligned bases of each of g's sequences, from its .hgn. */
void runsBySequence(const Genome *g, vector<uint64_t> &runs, vector<uint64_t> &bases) {
    lmap_file *f = g->hgn();
    uint32_t n = g->numSequences();
    runs.assign(n, 0);
    bases.assign(n, 0);
    for (uint32_t i = 0; i < n; i++) {
        lmap_cursor *c = lmap_cursor_open(f, 0, i, nullptr, 0, 0);
        if (!c) throw runtime_error("cannot read " + g->name() + "'s genome file");
        lmap_hit h;
        int rc = lmap_cursor_seek(c, 0, g->sequenceLength(i));
        while (rc == 0 && (rc = lmap_cursor_next(c, &h)) == 1) {
            rc = 0;
            runs[i]++;
            bases[i] += (uint64_t)h.len;
        }
        lmap_cursor_close(c);
        if (rc != 0) throw runtime_error("cannot read " + g->name() + "'s genome file");
    }
}

void summary(const Alignment &aln) {
    cout << "\nhal format 3\n" << aln.manifest().newick << "\n\n";
    cout << "GenomeName, NumChildren, Length, NumSequences, NumRuns, GenomeFileBytes, DnaBytes\n";
    for (const Genome *g : preorder(aln)) {
        const GenomeEntry *e = aln.manifest().find(g->name());
        lmap_file_stats st;
        lmap_get_stats(g->hgn(), &st);
        cout << g->name() << ", " << g->children().size() << ", " << e->length << ", " << e->numSequences << ", " << st.runs
             << ", " << e->hgnBytes << ", " << e->dnaBytes << "\n";
    }
    cout << "\n";
}

/* halStats' getGenomesInSpanningTree and printBranchPath */
void span(const Alignment &aln, const vector<string> &names, bool keepRoot) {
    vector<const Genome *> in;
    for (const string &n : names) in.push_back(need(aln, n));
    if (in.empty()) throw CliError("no genomes given");
    const Genome *lca = in[0];
    for (const Genome *g : in) lca = lowestCommonAncestor(lca, g);
    set<const Genome *> tree{lca};
    for (const Genome *g : in)
        for (const Genome *x = g; x != lca; x = x->parent()) tree.insert(x);
    vector<const Genome *> out;
    if (names.size() == 2) {
        /* the path from the first to the second */
        vector<const Genome *> up, down;
        for (const Genome *x = in[0]; x != lca; x = x->parent()) up.push_back(x);
        for (const Genome *x = in[1]; x != lca; x = x->parent()) down.push_back(x);
        out = up;
        out.push_back(lca);
        out.insert(out.end(), down.rbegin(), down.rend());
    } else {
        /* halStats prints a set of genome pointers, whose order is its memory layout;
         * tree preorder here */
        for (const Genome *g : preorder(aln))
            if (tree.count(g)) out.push_back(g);
    }
    for (const Genome *g : out)
        if (keepRoot || (g->parent() && tree.count(g->parent()))) cout << g->name() << " ";
    cout << endl;
}

unique_ptr<tb_file, void (*)(tb_file *)> openDna(const Alignment &aln, const Genome *g) {
    return unique_ptr<tb_file, void (*)(tb_file *)>(aln.openDna(g), tb_close);
}

/* halStats --baseComp: every step-th base of the concatenated genome */
void baseComp(const Alignment &aln, const string &arg) {
    vector<string> t = split(arg, ',');
    int64_t step = 0;
    if (t.size() == 2) step = atoll(t[1].c_str());
    if (step <= 0) throw CliError("Invalid value for --baseComp: " + arg + ".  Must be of format genomeName,step");
    const Genome *g = need(aln, t[0]);
    auto dna = openDna(aln, g);
    int64_t len = g->sequenceStart(g->numSequences() - 1) + g->sequenceLength(g->numSequences() - 1);
    if (step >= len) step = len - 1;
    uint64_t n[4] = {0, 0, 0, 0};
    uint32_t seq = 0;
    for (int64_t i = 0; i < len; i += step) {
        while (i >= g->sequenceStart(seq) + g->sequenceLength(seq)) seq++;
        char c;
        int64_t ti = tb_find(dna.get(), g->sequenceName(seq)), p = i - g->sequenceStart(seq);
        if (ti < 0 || tb_fetch(dna.get(), (uint32_t)ti, p, p + 1, &c, TB_NOMASK) != 0)
            throw runtime_error("cannot read DNA of " + g->name());
        const char *k = strchr("ACGT", c);
        if (k && c) n[k - "ACGT"]++;
    }
    double total = (double)(n[0] + n[1] + n[2] + n[3]);
    cout << n[0] / total << '\t' << n[1] / total << '\t' << n[2] / total << '\t' << n[3] / total << '\n';
}

/* halStats --assembly, figure for figure (including its integer mean and its median) */
void assembly(const Alignment &aln) {
    cout << "Genome, TotalLength, NumContigs, MinContigLength, MaxContigLength, MeanContigLength, MedianContigLength, N50, L50"
         << endl;
    for (const Genome *g : preorder(aln)) {
        multiset<int64_t> lengths;
        int64_t total = 0;
        for (uint32_t i = 0; i < g->numSequences(); i++) {
            lengths.insert(g->sequenceLength(i));
            total += g->sequenceLength(i);
        }
        int64_t minL = 0, maxL = 0, avg = 0, med = 0;
        if (!lengths.empty()) {
            minL = *lengths.begin();
            maxL = *lengths.rbegin();
            avg = (int64_t)((double)total / (double)lengths.size());
        }
        int64_t running = 0, count = 0, l50 = 0, n50 = 0;
        for (auto j = lengths.rbegin(); j != lengths.rend(); ++j) {
            running += *j;
            ++count;
            if (l50 == 0 && running >= total / 2) {
                n50 = *j;
                l50 = count;
            }
            if (count == (int64_t)(lengths.size() / 2)) med = *j;
        }
        cout << g->name() << ", " << total << ", " << lengths.size() << ", " << minL << ", " << maxL << ", " << avg << ", "
             << med << ", " << n50 << ", " << l50 << endl;
    }
}
}

int main_stats(int argc, char **argv) {
    return runCommand("stats", [&]() {
        enum {
            GENOMES = 1000, TREE, ROOT, CHILDREN, PARENT, BRANCH_LENGTH, BRANCHES, SPAN, SPAN_ROOT, SEQUENCES, SEQUENCE_STATS,
            CHROM_SIZES, BED_SEQUENCES, RUNS, BASE_COMP, GENOME_META, META, ASSEMBLY, NUM_SEGMENTS, TOP_SEGMENTS,
            BOTTOM_SEGMENTS, PERCENT_ID, COVERAGE, ALL_COVERAGE
        };
        static const option longOpts[] = {{"genomes", no_argument, nullptr, GENOMES},
                                          {"tree", no_argument, nullptr, TREE},
                                          {"root", no_argument, nullptr, ROOT},
                                          {"children", required_argument, nullptr, CHILDREN},
                                          {"parent", required_argument, nullptr, PARENT},
                                          {"branchLength", required_argument, nullptr, BRANCH_LENGTH},
                                          {"branches", no_argument, nullptr, BRANCHES},
                                          {"span", required_argument, nullptr, SPAN},
                                          {"spanRoot", required_argument, nullptr, SPAN_ROOT},
                                          {"sequences", required_argument, nullptr, SEQUENCES},
                                          {"sequenceStats", required_argument, nullptr, SEQUENCE_STATS},
                                          {"chromSizes", required_argument, nullptr, CHROM_SIZES},
                                          {"bedSequences", required_argument, nullptr, BED_SEQUENCES},
                                          {"runs", required_argument, nullptr, RUNS},
                                          {"baseComp", required_argument, nullptr, BASE_COMP},
                                          {"genomeMetaData", required_argument, nullptr, GENOME_META},
                                          {"metaData", no_argument, nullptr, META},
                                          {"assembly", no_argument, nullptr, ASSEMBLY},
                                          {"numSegments", required_argument, nullptr, NUM_SEGMENTS},
                                          {"topSegments", required_argument, nullptr, TOP_SEGMENTS},
                                          {"bottomSegments", required_argument, nullptr, BOTTOM_SEGMENTS},
                                          {"percentID", required_argument, nullptr, PERCENT_ID},
                                          {"coverage", required_argument, nullptr, COVERAGE},
                                          {"allCoverage", no_argument, nullptr, ALL_COVERAGE},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int what = 0, c;
        string arg;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            if (c == 'h') {
                cout << STATS_USAGE;
                return 0;
            }
            if (c < GENOMES) throw CliError("bad option");
            if (c == NUM_SEGMENTS || c == TOP_SEGMENTS || c == BOTTOM_SEGMENTS)
                throw CliError("format 3 has no segments: see --runs, and the runs in the summary and --sequenceStats");
            if (c == PERCENT_ID || c == COVERAGE || c == ALL_COVERAGE)
                throw CliError("--percentID and --coverage are not implemented for format 3 yet");
            if (what) throw CliError("the options are exclusive: give at most one");
            what = c;
            if (optarg) arg = optarg;
        }
        if (argc - optind != 1) throw CliError("expected 1 argument, got " + to_string(argc - optind));
        Alignment aln(argv[optind]);
        switch (what) {
        case 0: summary(aln); break;
        case GENOMES: {
            vector<const Genome *> gs = preorder(aln);
            for (size_t i = 0; i < gs.size(); i++) cout << gs[i]->name() << (i + 1 < gs.size() ? " " : "");
            cout << endl;
            break;
        }
        case TREE: cout << aln.manifest().newick << endl; break;
        case ROOT: cout << aln.root()->name() << endl; break;
        case CHILDREN: {
            const vector<const Genome *> &ch = need(aln, arg)->children();
            for (size_t i = 0; i < ch.size(); i++) cout << ch[i]->name() << (i + 1 < ch.size() ? " " : "");
            cout << endl;
            break;
        }
        case PARENT: {
            const Genome *g = need(aln, arg);
            if (g->parent()) cout << g->parent()->name() << endl;
            break;
        }
        case BRANCH_LENGTH: {
            const Genome *g = need(aln, arg);
            if (g->parent()) cout << aln.manifest().find(g->name())->branchLength << endl;
            break;
        }
        case BRANCHES: {
            bool first = true;
            for (const Genome *g : preorder(aln)) {
                if (!g->parent()) continue;
                cout << (first ? "" : " ") << g->name();
                first = false;
            }
            cout << endl;
            break;
        }
        case SPAN: span(aln, split(arg, ','), false); break;
        case SPAN_ROOT: span(aln, split(arg, ','), true); break;
        case SEQUENCES: {
            const Genome *g = need(aln, arg);
            for (uint32_t i = 0; i < g->numSequences(); i++) cout << (i ? "," : "") << g->sequenceName(i);
            cout << endl;
            break;
        }
        case SEQUENCE_STATS: {
            const Genome *g = need(aln, arg);
            vector<uint64_t> runs, bases;
            runsBySequence(g, runs, bases);
            if (g->numSequences() > 0) cout << "SequenceName, Length, NumRuns, AlignedToParent" << endl;
            for (uint32_t i = 0; i < g->numSequences(); i++)
                cout << g->sequenceName(i) << ", " << g->sequenceLength(i) << ", " << runs[i] << ", " << bases[i] << "\n";
            cout << endl;
            break;
        }
        case CHROM_SIZES: {
            const Genome *g = need(aln, arg);
            for (uint32_t i = 0; i < g->numSequences(); i++) cout << g->sequenceName(i) << '\t' << g->sequenceLength(i) << '\n';
            break;
        }
        case BED_SEQUENCES: {
            const Genome *g = need(aln, arg);
            for (uint32_t i = 0; i < g->numSequences(); i++)
                cout << g->sequenceName(i) << "\t" << 0 << "\t" << g->sequenceLength(i) << "\n";
            break;
        }
        case RUNS: {
            const Genome *g = need(aln, arg);
            lmap_file *f = g->hgn();
            for (uint32_t i = 0; i < g->numSequences(); i++) {
                lmap_cursor *cur = lmap_cursor_open(f, 0, i, nullptr, 0, 0);
                if (!cur) throw runtime_error("cannot read " + g->name() + "'s genome file");
                lmap_hit h;
                int rc = lmap_cursor_seek(cur, 0, g->sequenceLength(i));
                while (rc == 0 && (rc = lmap_cursor_next(cur, &h)) == 1) {
                    rc = 0;
                    cout << g->sequenceName(i) << '\t' << h.a << '\t' << h.a + h.len << '\t'
                         << lmap_member_at(f, 1, h.b_member)->name << ':' << h.b << "\t0\t" << (h.strand ? '-' : '+')
                         << '\n';
                }
                lmap_cursor_close(cur);
                if (rc != 0) throw runtime_error("cannot read " + g->name() + "'s genome file");
            }
            break;
        }
        case BASE_COMP: baseComp(aln, arg); break;
        case GENOME_META: {
            lmap_file *f = need(aln, arg)->hgn();
            const string prefix = "hal.meta.";
            for (uint32_t i = 0; i < lmap_meta_count(f); i++) {
                const char *k, *v;
                if (lmap_meta_at(f, i, &k, &v) == 0 && !strncmp(k, prefix.c_str(), prefix.size()))
                    cout << k + prefix.size() << '\t' << v << '\n';
            }
            break;
        }
        case META:
            for (const auto &kv : aln.manifest().meta) cout << kv.first << '\t' << kv.second << '\n';
            break;
        case ASSEMBLY: assembly(aln); break;
        }
        cout.flush();
        if (!cout) throw runtime_error("write failed");
        return 0;
    });
}
}
