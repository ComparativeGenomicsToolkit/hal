/*
 * hal 4d: fourfold degenerate sites of annotated coding sequence, as universal columns,
 * for training a neutral model (phyloFit) on the alignment's own columns.
 *
 * Each annotation (GFF3 or GTF, gzipped or not: its CDS records) gives codons of its
 * genome.  A codon whose first two bases are a fourfold degenerate prefix there is
 * followed into the alignment: each of its three bases' columns -- the base's reference
 * column, duplicates kept, which is the universal column it lies in -- and it is kept as
 * msa_view --4d keeps one: every leaf with a base in any of the three columns has a base in
 * all three, the first two positions agree in all of them, and (as phyloFit has one row
 * per genome) no genome, leaf or ancestor, has two copies in any of the three.  Its site is
 * the third position's column.  A column is identified by its top base, so one chosen by
 * several annotations counts once; one their codons (any, not only fourfold degenerate ones)
 * give different roles (codon position or strand) is dropped.  A codon whose site's column an
 * earlier annotation tested, or (--sample) that hashes above the sample as its annotation
 * starts, is not tested: the work after the first annotation is the new columns' alone.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Liftover.h"
#include "hal3Walk.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <zlib.h>

using namespace std;

namespace hal3 {

namespace {

const char *FOURD_USAGE =
    "usage: hal 4d [options] <alignment> <out.ss> <genome>=<annotation> [<genome>=<annotation> ...]\n"
    "Fourfold degenerate sites of annotated coding sequence, in universal columns, as a phast\n"
    "SS file (tuple size 1, one row per leaf, phast names: a genome's name up to its first\n"
    "dot) for phyloFit.  Each annotation is a GFF3 or GTF file (gzipped or not) of a leaf\n"
    "genome; its CDS records, grouped by transcript (Parent=, or transcript_id) and framed by\n"
    "their phases, give codons.  A codon whose first two bases are a fourfold degenerate\n"
    "prefix is kept, as msa_view --4d keeps one, if every leaf with a base in any of its three\n"
    "columns has one in all three and the first two positions agree in all of them, with\n"
    "duplicated genomes handled by --dupes.  Its site is the third position's column (or a\n"
    "single-copy part of it), in the codon's orientation.  Columns are hal maf's reference\n"
    "columns with duplicates kept, which are the universal columns, identified by their top\n"
    "base (a part by its root): one chosen by several annotations (or transcripts) counts\n"
    "once: a later annotation's codon whose site is a column already tested is not tested again\n"
    "(the test is the columns', not the codon's).\n"
    "Skipped: transcripts whose phases disagree, positions two transcripts of an annotation put\n"
    "in different frames, and columns the annotations' codons (any, not only fourfold\n"
    "degenerate ones) give different roles (codon position or strand).\n"
    "options:\n"
    "  --dupes split|mask|drop  a genome with two copies in a codon's columns.  split (the\n"
    "                       default): the third position's column is cut where a base has two\n"
    "                       copies of a genome below it, into single-copy parts, each a site\n"
    "                       (its leaves' first two bases their own copies', at the neighbouring\n"
    "                       positions: codons across an intron are skipped); a part needs two\n"
    "                       leaves, and the column's site is its largest part that passes (the\n"
    "                       most leaves; a tie to the lower hash of its root), or with --parts all\n"
    "                       every part that does.  mask: the genome's leaves (an ancestor's: all\n"
    "                       below it) left out, the codon dropped if the annotated genome is\n"
    "                       among them.  drop: the codon dropped\n"
    "  --dupStats FILE      per genome: codons tested where it has two copies\n"
    "  --sample N           keep N sites: those whose columns hash lowest (deterministic, and a\n"
    "                       smaller sample is part of a larger one).  Only the 2N lowest are\n"
    "                       held; a codon whose site hashes above them as its annotation starts\n"
    "                       is not tested (except with --parts all)\n"
    "  --sites FILE         each site held (with --sample N, the 2N that hash lowest): its\n"
    "                       column's top base (or its part's root), the codon's strand against\n"
    "                       the column's top base, how many annotations chose it (the first that\n"
    "                       kept it, and later ones with a codon whose site is in it: not tested\n"
    "                       again, except with --parts all), the first, and what became of it\n"
    "  --stats FILE         per annotation: transcripts, codons and why codons were dropped (TSV)\n"
    "  --excludeGenomes A,B leaves left out (no row, and not walked into)\n"
    "  --threads N          N threads (default 1); the output is the same\n";

/* a codon: its bases' positions in transcription order, on sequence seq */
struct Codon {
    uint32_t seq;
    int64_t p[3];
    bool minus;
    bool operator<(const Codon &o) const {
        return tie(seq, p[0], p[1], p[2], minus) < tie(o.seq, o.p[0], o.p[1], o.p[2], o.minus);
    }
    bool operator==(const Codon &o) const { return seq == o.seq && p[0] == o.p[0] && p[1] == o.p[1] && p[2] == o.p[2] && minus == o.minus; }
};

struct AnnStats {
    uint64_t cdsLines = 0, transcripts = 0, unknownSeq = 0, phaseBad = 0, codons = 0, frameConflict = 0;
    uint64_t anchorN = 0, notPrefix = 0, candidates = 0, dupLeaf = 0, dupAnc = 0, dupAnchor = 0, masked = 0, gap = 0, differs = 0,
             passed = 0, alone = 0, spliced = 0, parts = 0, partsAlone = 0, widest = 0, columns = 0, held = 0, tried = 0, above = 0;
    double liftSecs = 0, sweepSecs = 0, testSecs = 0;   /* thread time in each step */
    void add(const AnnStats &o) {
        liftSecs += o.liftSecs;
        sweepSecs += o.sweepSecs;
        testSecs += o.testSecs;
        anchorN += o.anchorN;
        notPrefix += o.notPrefix;
        candidates += o.candidates;
        held += o.held;
        tried += o.tried;
        above += o.above;
        dupLeaf += o.dupLeaf;
        dupAnc += o.dupAnc;
        dupAnchor += o.dupAnchor;
        masked += o.masked;
        spliced += o.spliced;
        columns += o.columns;
        widest = max(widest, o.widest);
        parts += o.parts;
        partsAlone += o.partsAlone;
        gap += o.gap;
        differs += o.differs;
        passed += o.passed;
        alone += o.alone;
    }
};

/* a column, by its top base */
struct ColKey {
    uint32_t genome, seq;
    int64_t pos;
    bool operator==(const ColKey &o) const { return genome == o.genome && seq == o.seq && pos == o.pos; }
    bool operator<(const ColKey &o) const { return tie(genome, seq, pos) < tie(o.genome, o.seq, o.pos); }
};

struct ColKeyHash {
    size_t operator()(const ColKey &k) const {
        uint64_t h = (uint64_t)k.genome * 0x9E3779B97F4A7C15ULL ^ (uint64_t)k.seq * 0xC2B2AE3D27D4EB4FULL ^ (uint64_t)k.pos * 0x165667B19E3779F9ULL;
        return (size_t)(h ^ (h >> 29));
    }
};

/* a codon kept: its third position's column and orientation, its site, and the site's bases
 * (by leaf) */
struct Kept {
    ColKey col;      /* the third position's column, by its top base */
    bool minus;      /* the codon's strand, against that top base (shared by all annotations) */
    ColKey site;     /* the site: its column (drop, mask, split), or its single-copy part's root (--parts all) */
    uint64_t hash;   /* the site's, for sampling */
    string tuple;    /* empty: the site was held as the annotation started (it counts the annotation, untested) */
};

/* a role a codon gives a column: bit 2 x position + (the codon on the minus strand of the
 * column's top base) */
struct Role {
    ColKey col;
    uint64_t hash;   /* the column's, as a site's */
    uint8_t bit;
};

bool gzLine(gzFile f, string &line) {
    line.clear();
    char buf[1 << 16];
    while (gzgets(f, buf, sizeof buf)) {
        line += buf;
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
            if (!line.empty() && line.back() == '\r') line.pop_back();
            return true;
        }
    }
    return !line.empty();
}

/* attribute key's values: GFF3 key=a,b or GTF key "a" */
vector<string> attribute(const string &attrs, const string &key) {
    vector<string> out;
    size_t a = 0;
    while (a < attrs.size()) {
        size_t b = attrs.find(';', a);
        if (b == string::npos) b = attrs.size();
        size_t s = a;
        while (s < b && attrs[s] == ' ') s++;
        if (attrs.compare(s, key.size(), key) == 0 && s + key.size() < b) {
            char c = attrs[s + key.size()];
            if (c == '=') {
                string v = attrs.substr(s + key.size() + 1, b - s - key.size() - 1);
                stringstream ss(v);
                string x;
                while (getline(ss, x, ','))
                    if (!x.empty()) out.push_back(x);
            } else if (c == ' ') {
                size_t q0 = attrs.find('"', s + key.size()), q1 = q0 == string::npos ? q0 : attrs.find('"', q0 + 1);
                if (q0 != string::npos && q1 != string::npos && q1 < b) out.push_back(attrs.substr(q0 + 1, q1 - q0 - 1));
            }
        }
        a = b + 1;
    }
    return out;
}

/* the codons of an annotation's CDS, on genome g: unique, without frame conflicts */
vector<Codon> readCodons(const string &path, const Genome *g, AnnStats &st) {
    gzFile f = gzopen(path.c_str(), "r");
    if (!f) throw runtime_error("cannot read " + path);
    struct Seg {
        int64_t seq, start, end;
        char strand;
        int phase;
    };
    map<string, vector<Seg>> tx;
    string line;
    vector<string> fld;
    while (gzLine(f, line)) {
        if (line.empty() || line[0] == '#') {
            if (line.compare(0, 7, "##FASTA") == 0) break;
            continue;
        }
        fld.clear();
        size_t a = 0;
        for (int i = 0; i < 9; i++) {
            size_t b = i < 8 ? line.find('\t', a) : line.size();
            if (b == string::npos) break;
            fld.push_back(line.substr(a, b - a));
            a = b + 1;
        }
        if (fld.size() < 9 || fld[2] != "CDS") continue;
        st.cdsLines++;
        Seg s{g->findSequence(fld[0]), atoll(fld[3].c_str()) - 1, atoll(fld[4].c_str()), fld[6].empty() ? '.' : fld[6][0],
              fld[7] == "." ? 0 : atoi(fld[7].c_str())};
        vector<string> ids = attribute(fld[8], "Parent");
        if (ids.empty()) ids = attribute(fld[8], "transcript_id");
        if (ids.empty()) ids = attribute(fld[8], "ID");
        for (const string &id : ids) tx[id].push_back(s);
    }
    gzclose(f);

    vector<Codon> codons;
    for (auto &t : tx) {
        vector<Seg> &v = t.second;
        st.transcripts++;
        bool unknown = false, mixed = false;
        for (const Seg &s : v) {
            unknown = unknown || s.seq < 0;
            mixed = mixed || s.seq != v[0].seq || s.strand != v[0].strand || (s.strand != '+' && s.strand != '-') || s.end <= s.start ||
                    s.phase < 0 || s.phase > 2;
        }
        if (unknown) {
            st.unknownSeq++;
            continue;
        }
        if (mixed) {
            st.phaseBad++;
            continue;
        }
        const bool minus = v[0].strand == '-';
        sort(v.begin(), v.end(), [&](const Seg &x, const Seg &y) { return minus ? x.end > y.end : x.start < y.start; });
        /* a segment's phase: what of it ends the previous codon */
        int64_t cum = 0;
        const int64_t phase0 = v[0].phase;
        bool ok = true;
        for (size_t k = 0; k < v.size() && ok; k++) {
            if (k > 0 && ((phase0 - cum) % 3 + 3) % 3 != v[k].phase) ok = false;
            if (k > 0 && (minus ? v[k].end > v[k - 1].start : v[k].start < v[k - 1].end)) ok = false;   /* overlapping */
            cum += v[k].end - v[k].start;
        }
        if (!ok) {
            st.phaseBad++;
            continue;
        }
        /* codon positions along the concatenated CDS */
        auto at = [&](int64_t t) {
            for (const Seg &s : v) {
                const int64_t len = s.end - s.start;
                if (t < len) return minus ? s.end - 1 - t : s.start + t;
                t -= len;
            }
            return (int64_t)-1;
        };
        for (int64_t t = phase0; t + 2 < cum; t += 3) codons.push_back(Codon{(uint32_t)v[0].seq, {at(t), at(t + 1), at(t + 2)}, minus});
    }
    sort(codons.begin(), codons.end());
    codons.erase(unique(codons.begin(), codons.end()), codons.end());
    st.codons = codons.size();
    /* positions the annotation puts in two frames (or on two strands) */
    vector<pair<uint64_t, uint8_t>> roles;
    roles.reserve(3 * codons.size());
    for (const Codon &c : codons)
        for (int i = 0; i < 3; i++) roles.push_back(make_pair((uint64_t)c.seq << 40 | (uint64_t)c.p[i], (uint8_t)(i | (c.minus ? 4 : 0))));
    sort(roles.begin(), roles.end());
    vector<uint64_t> conflicted;
    for (size_t i = 1; i < roles.size(); i++)
        if (roles[i].first == roles[i - 1].first && roles[i].second != roles[i - 1].second &&
            (conflicted.empty() || conflicted.back() != roles[i].first))
            conflicted.push_back(roles[i].first);
    vector<Codon> out;
    for (const Codon &c : codons) {
        bool bad = false;
        for (int i = 0; i < 3 && !bad; i++) bad = binary_search(conflicted.begin(), conflicted.end(), (uint64_t)c.seq << 40 | (uint64_t)c.p[i]);
        if (bad) st.frameConflict++;
        else out.push_back(c);
    }
    return out;
}

char complementBase(char c) {
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: return 'N';
    }
}

char upperBase(char c) {
    c = (char)toupper((unsigned char)c);
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' ? c : 'N';
}

/* a fourfold degenerate first two bases: Ser TCN, Thr ACN, Leu CTN, Pro CCN, Arg CGN, Val GTN,
 * Ala GCN, Gly GGN */
bool fourfold(char a, char b) {
    return ((a == 'A' || a == 'T') && b == 'C') || ((a == 'C' || a == 'G') && b != 'A' && b != 'N');
}

/* one unit of work: an annotation's codons on one sequence, close together */
struct Unit {
    size_t first, last;   /* codons [first, last) */
};

/* a base of a column: its genome, its parent (an entry at the same position), where it is,
 * its orientation against the anchor's, and (a leaf's) its base in the anchor's orientation */
struct Entry {
    uint32_t genome;
    int32_t parent;
    uint32_t seq;
    int64_t t;
    bool rev;
    char base;
};

/* a column's entries */
struct Col {
    const Entry *b, *e;
    size_t size() const { return (size_t)(e - b); }
    bool empty() const { return b == e; }
    const Entry &operator[](size_t i) const { return b[i]; }
    const Entry *begin() const { return b; }
    const Entry *end() const { return e; }
};

/* a lift piece's positions [i0, i1) (indices into the unit's positions); its entries' slots
 * (each counted from its column's first) from slots[start] */
struct Span {
    size_t i0, i1, start;
};

/* what a thread keeps between units: its walk (and cursors) for the current anchor, DNA, and
 * the units' buffers */
struct Worker {
    const Genome *anchor = nullptr;
    unique_ptr<TreeWalk> walk;
    unique_ptr<DnaReader> dna;
    vector<Entry> ent;                  /* every column's entries, column after column */
    vector<size_t> colStart, fill;      /* per position: its column's first entry; entries placed */
    vector<int32_t> slots, orderOf;     /* orderOf: genome index -> its place in the unit's order, or -1 */
    vector<Span> spans;
    vector<size_t> spanStart;           /* per genome in order: its pieces' first span */
    vector<uint64_t> kidKeys;
    vector<int32_t> root;
    vector<const Entry *> at[2];
    vector<pair<int32_t, size_t>> parts;
};

struct UnitResult {
    vector<Kept> kept;
    vector<Role> roles;
    vector<pair<ColKey, uint64_t>> tested;   /* the site columns tested (not --parts all), and their hashes */
    AnnStats st;
    map<uint32_t, uint64_t> dupBy;   /* genome -> candidate codons where it has two copies */
};

struct Context {
    const Alignment *aln;
    WalkOptions walkOpts;
    map<const Genome *, uint32_t> index;   /* every genome */
    vector<int32_t> leafOf;                /* genome index -> leaf row, -1 */
    vector<pair<int32_t, int32_t>> rows;   /* genome index -> its leaves' rows [lo, hi) */
    vector<int> depth;                     /* genome index -> its depth (the root 0) */
    size_t leaves;
    enum Dupes { DROP, MASK, SPLIT } dupes = SPLIT;
    bool allParts = false;                 /* split: every part a site (else a column's largest) */
    vector<char> excluded;                 /* genome index -> left out */
    vector<const Genome *> genomes;        /* by index */
    /* as the annotation starts: the sample's cutoff (--sample: a site hashing higher can no
     * longer be held), the sites held, and the site columns tested (not --parts all: a codon
     * whose site is one is not tested again; the test is the columns', not the codon's) */
    uint64_t cutoff = UINT64_MAX;
    unordered_set<ColKey, ColKeyHash> held;
    unordered_map<ColKey, uint64_t, ColKeyHash> tried;   /* and its hash */
};

uint64_t fnv(const string &s, uint64_t h = 1469598103934665603ULL);

uint64_t siteHash(const Context &cx, const ColKey &k) {
    const Genome *g = cx.genomes[k.genome];
    return fnv(g->name() + '\0' + g->sequenceName(k.seq) + '\0' + to_string(k.pos));
}

UnitResult runUnit(const Context &cx, Worker &wk, const Genome *anchor, const vector<Codon> &codons, const Unit &u) {
    UnitResult r;
    if (wk.anchor != anchor) {
        wk.walk.reset(new TreeWalk(*cx.aln, anchor, cx.walkOpts));
        if (!wk.dna) wk.dna.reset(new DnaReader(*cx.aln));
        wk.anchor = anchor;
    }
    const uint32_t seq = codons[u.first].seq;
    auto intervals = [](const vector<int64_t> &v) {   /* sorted, unique positions */
        vector<pair<int64_t, int64_t>> iv;
        for (int64_t x : v) {
            if (!iv.empty() && iv.back().second == x) iv.back().second = x + 1;
            else iv.push_back(make_pair(x, x + 1));
        }
        return iv;
    };
    /* every codon's positions, and their columns' top bases (the way up alone) */
    vector<int64_t> all;
    for (size_t i = u.first; i < u.last; i++)
        for (int k = 0; k < 3; k++) all.push_back(codons[i].p[k]);
    sort(all.begin(), all.end());
    all.erase(unique(all.begin(), all.end()), all.end());
    auto t0 = chrono::steady_clock::now();
    PieceMap tops;
    wk.walk->liftTops(seq, intervals(all), tops);
    vector<ColKey> allTop(all.size());
    vector<char> allRev(all.size(), 2);
    for (auto &tp : tops) {
        const uint32_t gi = cx.index.at(tp.first);
        for (const LiftPiece &p : tp.second)
            for (size_t i = (size_t)(lower_bound(all.begin(), all.end(), p.srcStart) - all.begin()); i < all.size() && all[i] < p.srcStart + p.length; i++) {
                const int64_t d = all[i] - p.srcStart;
                allTop[i] = ColKey{gi, p.tgtSeq, p.reversed ? p.tgtStart + p.length - 1 - d : p.tgtStart + d};
                allRev[i] = p.reversed ? 1 : 0;
            }
    }
    vector<uint64_t> allHash(all.size());
    for (size_t i = 0; i < all.size(); i++) {
        if (allRev[i] == 2) throw runtime_error("internal error: hal 4d: a position with no top");
        allHash[i] = siteHash(cx, allTop[i]);
    }
    auto allAt = [&](int64_t x) { return (size_t)(lower_bound(all.begin(), all.end(), x) - all.begin()); };
    /* the roles every codon gives its columns; with --sample, of the columns that can still
     * be sites (with --parts all a site is a part, and its column's hash says nothing) */
    for (size_t i = u.first; i < u.last; i++) {
        const Codon &c = codons[i];
        for (int k = 0; k < 3; k++) {
            const size_t j = allAt(c.p[k]);
            if (cx.allParts || allHash[j] <= cx.cutoff)
                r.roles.push_back(Role{allTop[j], allHash[j], (uint8_t)(1 << (2 * k + (c.minus != (allRev[j] != 0) ? 1 : 0)))});
        }
    }

    /* the anchor's bases over the unit: its codons' prefixes */
    string adna;
    const int64_t lo = all.front();
    wk.dna->fetch(anchor, seq, lo, all.back() + 1, adna, true);   /* upper-cased anyway: no soft-masking */
    auto coding = [](const Codon &c, char forward) { return c.minus ? complementBase(upperBase(forward)) : upperBase(forward); };
    vector<const Codon *> cand;
    vector<int64_t> pos;
    for (size_t i = u.first; i < u.last; i++) {
        const Codon &c = codons[i];
        char a = coding(c, adna[(size_t)(c.p[0] - lo)]), b = coding(c, adna[(size_t)(c.p[1] - lo)]), t = coding(c, adna[(size_t)(c.p[2] - lo)]);
        if (a == 'N' || b == 'N' || t == 'N') {
            r.st.anchorN++;
            continue;
        }
        if (!fourfold(a, b)) {
            r.st.notPrefix++;
            continue;
        }
        r.st.candidates++;
        /* split: a part's leaves' first two bases are their own copies', at the neighbouring
         * positions of that copy, so a codon across an intron cannot be followed */
        if (cx.dupes == Context::SPLIT && !((c.p[1] == c.p[0] + 1 && c.p[2] == c.p[1] + 1) || (c.p[1] == c.p[0] - 1 && c.p[2] == c.p[1] - 1))) {
            r.st.spliced++;
            continue;
        }
        if (!cx.allParts) {   /* its site is its third position's column */
            const size_t j = allAt(c.p[2]);
            if (cx.held.count(allTop[j])) {   /* held: the annotation counts, untested */
                r.st.held++;
                Kept k;
                k.site = allTop[j];
                k.hash = allHash[j];
                r.kept.push_back(move(k));
                continue;
            }
            if (cx.tried.count(allTop[j])) {   /* tested before, and not held */
                r.st.tried++;
                continue;
            }
            if (allHash[j] > cx.cutoff) {
                r.st.above++;
                continue;
            }
        }
        cand.push_back(&c);
        for (int k = 0; k < 3; k++) pos.push_back(c.p[k]);
    }
    auto t1 = chrono::steady_clock::now();
    r.st.liftSecs += chrono::duration<double>(t1 - t0).count();
    if (cand.empty()) return r;
    sort(pos.begin(), pos.end());
    pos.erase(unique(pos.begin(), pos.end()), pos.end());
    PieceMap pieces;
    map<const Genome *, vector<int32_t>> from;
    wk.walk->liftParented(seq, intervals(pos), pieces, from);
    const auto tl = chrono::steady_clock::now();
    r.st.liftSecs += chrono::duration<double>(tl - t1).count();
    t1 = tl;

    /* each position's column as entries, genomes in preorder (a parent before its
     * children), in wk.ent column after column: counted first, then placed */
    const size_t n = pos.size();
    const uint32_t ai = cx.index.at(anchor);
    vector<pair<uint32_t, const Genome *>> order;
    for (auto &gp : pieces) {
        auto gi = cx.index.find(gp.first);
        if (gi != cx.index.end() && !cx.excluded[gi->second] && !gp.second.empty()) order.push_back(make_pair(gi->second, gp.first));
    }
    sort(order.begin(), order.end());
    wk.orderOf.assign(cx.genomes.size(), -1);
    wk.spans.clear();
    wk.spanStart.clear();
    wk.colStart.assign(n + 1, 0);
    for (size_t o = 0; o < order.size(); o++) {
        wk.orderOf[order[o].first] = (int32_t)o;
        wk.spanStart.push_back(wk.spans.size());
        for (const LiftPiece &p : pieces.at(order[o].second)) {
            size_t i0 = (size_t)(lower_bound(pos.begin(), pos.end(), p.srcStart) - pos.begin()), i1 = i0;
            while (i1 < n && pos[i1] < p.srcStart + p.length) i1++;
            wk.spans.push_back(Span{i0, i1, 0});
            for (size_t i = i0; i < i1; i++) wk.colStart[i + 1]++;
        }
    }
    for (size_t i = 0; i < n; i++) wk.colStart[i + 1] += wk.colStart[i];
    wk.ent.resize(wk.colStart[n]);
    wk.fill.assign(n, 0);
    wk.slots.clear();
    auto col = [&](size_t i) { return Col{wk.ent.data() + wk.colStart[i], wk.ent.data() + wk.colStart[i + 1]}; };
    string dna;
    for (size_t o = 0; o < order.size(); o++) {
        const uint32_t gidx = order[o].first;
        const Genome *g = order[o].second;
        const int32_t leaf = cx.leafOf[gidx];
        const vector<LiftPiece> &ps = pieces.at(g);
        const vector<int32_t> &fs = from.at(g);
        const int32_t po = g->parent() ? wk.orderOf[cx.index.at(g->parent())] : -1;   /* the parent's place, if it has pieces */
        const size_t pn = po >= 0 ? wk.spanStart[(size_t)po + 1] - wk.spanStart[(size_t)po] : 0;
        for (size_t j = 0; j < ps.size(); j++) {
            const LiftPiece &p = ps[j];
            Span &sp = wk.spans[wk.spanStart[o] + j];
            sp.start = wk.slots.size();
            const size_t i0 = sp.i0, i1 = sp.i1;
            if (i0 == i1) continue;
            int64_t tlo = 0;
            if (leaf >= 0) {   /* the piece's bases at these positions, one fetch */
                int64_t a = p.reversed ? p.tgtStart + p.length - 1 - (pos[i1 - 1] - p.srcStart) : p.tgtStart + (pos[i0] - p.srcStart);
                int64_t b = p.reversed ? p.tgtStart + p.length - 1 - (pos[i0] - p.srcStart) : p.tgtStart + (pos[i1 - 1] - p.srcStart);
                wk.dna->fetch(g, p.tgtSeq, a, b + 1, dna, true);
                tlo = a;
            }
            const int32_t f = fs[j];
            if (f >= 0 && (size_t)f >= pn) throw runtime_error("internal error: hal 4d: a piece without its parent piece");
            const Span *q = f >= 0 ? &wk.spans[wk.spanStart[(size_t)po] + (size_t)f] : nullptr;
            for (size_t i = i0; i < i1; i++) {
                const int64_t d = pos[i] - p.srcStart, t = p.reversed ? p.tgtStart + p.length - 1 - d : p.tgtStart + d;
                int32_t parent = -1;
                if (q) {
                    if (i < q->i0 || i >= q->i1) throw runtime_error("internal error: hal 4d: a base outside its parent piece");
                    parent = wk.slots[q->start + i - q->i0];
                }
                char c = 0;
                if (leaf >= 0) {
                    c = upperBase(dna[(size_t)(t - tlo)]);
                    if (p.reversed) c = complementBase(c);
                }
                const size_t slot = wk.fill[i]++;
                wk.slots.push_back((int32_t)slot);
                wk.ent[wk.colStart[i] + slot] = Entry{gidx, parent, p.tgtSeq, t, p.reversed, c};
            }
        }
    }
    /* each position: its top (the shallowest base), and the genomes with two copies */
    vector<ColKey> top(n);
    vector<char> topRev(n, 0);
    vector<vector<uint32_t>> dupG(n);
    vector<uint32_t> seen(cx.leafOf.size(), UINT32_MAX);
    for (size_t i = 0; i < n; i++) {
        int best = INT32_MAX;
        for (const Entry &e : col(i)) {
            const int dep = (int)cx.depth[e.genome];
            if (dep < best) {
                best = dep;
                top[i] = ColKey{e.genome, e.seq, e.t};
                topRev[i] = e.rev;
            }
            if (seen[e.genome] == (uint32_t)i && (dupG[i].empty() || dupG[i].back() != e.genome)) dupG[i].push_back(e.genome);
            seen[e.genome] = (uint32_t)i;
        }
        if (col(i).empty()) throw runtime_error("internal error: hal 4d: a position with no column");
        const size_t j = allAt(pos[i]);
        if (!(top[i] == allTop[j]) || topRev[i] != allRev[j]) throw runtime_error("internal error: hal 4d: a column's top differs from the way up's");
        r.st.widest = max(r.st.widest, (uint64_t)col(i).size());
    }
    auto t2 = chrono::steady_clock::now();
    r.st.sweepSecs += chrono::duration<double>(t2 - t1).count();

    const int32_t ar = cx.leafOf[ai];
    for (const Codon *c : cand) {
        size_t ix[3];
        for (int k = 0; k < 3; k++) ix[k] = (size_t)(lower_bound(pos.begin(), pos.end(), c->p[k]) - pos.begin());
        auto oriented = [&](char x) { return c->minus ? complementBase(x) : x; };
        vector<uint32_t> dups;
        for (int k = 0; k < 3; k++) dups.insert(dups.end(), dupG[ix[k]].begin(), dupG[ix[k]].end());
        sort(dups.begin(), dups.end());
        dups.erase(unique(dups.begin(), dups.end()), dups.end());
        for (uint32_t d : dups) r.dupBy[d]++;

        if (cx.dupes == Context::SPLIT) {
            /* the third position's column cut where a base has two copies of a genome below
             * it: single-copy parts, each a column of its own */
            const Col e3 = col(ix[2]);
            vector<uint64_t> &kids = wk.kidKeys;   /* (parent, genome) of each entry with a parent */
            kids.clear();
            for (const Entry &e : e3)
                if (e.parent >= 0) kids.push_back((uint64_t)(uint32_t)e.parent << 32 | e.genome);
            sort(kids.begin(), kids.end());
            vector<int32_t> &root = wk.root;
            root.resize(e3.size());
            for (size_t i = 0; i < e3.size(); i++) {
                const Entry &e = e3[i];
                bool cut = e.parent < 0;
                if (!cut) {
                    const uint64_t key = (uint64_t)(uint32_t)e.parent << 32 | e.genome;
                    auto lo = lower_bound(kids.begin(), kids.end(), key);
                    cut = lo + 1 < kids.end() && lo[1] == key;   /* two copies of the genome below the parent */
                }
                root[i] = cut ? (int32_t)i : root[(size_t)e.parent];
            }
            auto byPlace = [](const Entry *a, const Entry *b) { return tie(a->genome, a->seq, a->t, a->rev) < tie(b->genome, b->seq, b->t, b->rev); };
            for (int k = 0; k < 2; k++) {
                wk.at[k].clear();
                for (const Entry &e : col(ix[k])) wk.at[k].push_back(&e);
                sort(wk.at[k].begin(), wk.at[k].end(), byPlace);
            }
            auto find = [&](int k, const Entry &want) -> const Entry * {
                auto it = lower_bound(wk.at[k].begin(), wk.at[k].end(), &want, byPlace);
                return it != wk.at[k].end() && !byPlace(&want, *it) ? *it : nullptr;
            };
            vector<pair<int32_t, size_t>> &parts = wk.parts;   /* (part root, leaf entry), by part */
            parts.clear();
            for (size_t i = 0; i < e3.size(); i++)
                if (cx.leafOf[e3[i].genome] >= 0) parts.push_back(make_pair(root[i], i));
            sort(parts.begin(), parts.end());
            size_t bestB = 0, bestE = 0;   /* the largest part that passes: parts[bestB, bestE) */
            uint64_t bestHash = 0;
            for (size_t pb = 0, pe = 0; pb < parts.size(); pb = pe) {
                pe = pb + 1;
                while (pe < parts.size() && parts[pe].first == parts[pb].first) pe++;
                const size_t size = pe - pb;
                if (size < 2) {
                    r.st.partsAlone++;
                    continue;
                }
                r.st.parts++;
                char p0 = 'N', p1 = 'N';
                bool missing = false, differs = false;
                for (size_t m = pb; m < pe; m++) {
                    const Entry &x = e3[parts[m].second];
                    const int64_t sg = x.rev ? -1 : 1;
                    const Entry *y[2] = {nullptr, nullptr};
                    for (int k = 0; k < 2; k++) {
                        Entry want = x;
                        want.t = x.t + (c->p[k] - c->p[2]) * sg;
                        y[k] = find(k, want);
                    }
                    if (!y[0] || !y[1]) {
                        missing = true;
                        break;
                    }
                    char b0 = oriented(y[0]->base), b1 = oriented(y[1]->base);
                    if (b0 != 'N') {
                        if (p0 == 'N') p0 = b0;
                        else if (b0 != p0) differs = true;
                    }
                    if (b1 != 'N') {
                        if (p1 == 'N') p1 = b1;
                        else if (b1 != p1) differs = true;
                    }
                }
                if (missing) {
                    r.st.gap++;
                    continue;
                }
                if (differs || !fourfold(p0, p1)) {
                    r.st.differs++;
                    continue;
                }
                r.st.passed++;
                const Entry &rt = e3[(size_t)parts[pb].first];
                const ColKey partRoot{rt.genome, rt.seq, rt.t};
                const uint64_t h = siteHash(cx, partRoot);
                if (!cx.allParts) {   /* the most leaves; a tie to the lower hash of its root */
                    if (bestE == bestB || size > bestE - bestB || (size == bestE - bestB && h < bestHash)) {
                        bestB = pb;
                        bestE = pe;
                        bestHash = h;
                    }
                    continue;
                }
                if (h > cx.cutoff) continue;   /* can no longer be held */
                Kept k;
                k.col = top[ix[2]];
                k.minus = c->minus != (topRev[ix[2]] != 0);
                k.site = partRoot;
                k.hash = h;
                k.tuple.assign(cx.leaves, '*');
                for (size_t m = pb; m < pe; m++) k.tuple[(size_t)cx.leafOf[e3[parts[m].second].genome]] = oriented(e3[parts[m].second].base);
                r.kept.push_back(move(k));
            }
            if (!cx.allParts) r.tested.push_back(make_pair(top[ix[2]], allHash[allAt(c->p[2])]));
            if (bestE > bestB) {   /* the column's site: its largest part */
                Kept k;
                k.col = k.site = top[ix[2]];
                k.minus = c->minus != (topRev[ix[2]] != 0);
                k.hash = allHash[allAt(c->p[2])];
                k.tuple.assign(cx.leaves, '*');
                for (size_t m = bestB; m < bestE; m++) k.tuple[(size_t)cx.leafOf[e3[parts[m].second].genome]] = oriented(e3[parts[m].second].base);
                r.kept.push_back(move(k));
                r.st.columns++;
            }
            continue;
        }

        /* drop, mask: the whole column, a duplicated genome's leaves (itself, or all below an
         * ancestor) having their bases from different copies */
        vector<pair<int32_t, int32_t>> masked;
        bool anchorDup = false, leafDup = false;
        for (uint32_t d : dups) {
            masked.push_back(cx.rows[d]);
            anchorDup = anchorDup || (ar >= cx.rows[d].first && ar < cx.rows[d].second);
            leafDup = leafDup || cx.leafOf[d] >= 0;
        }
        /* tested, but for the annotated genome's own copies: another's codon there is tested */
        if (!cx.allParts && (cx.dupes == Context::DROP || !anchorDup)) r.tested.push_back(make_pair(top[ix[2]], allHash[allAt(c->p[2])]));
        if (!dups.empty() && cx.dupes == Context::DROP) {
            (leafDup ? r.st.dupLeaf : r.st.dupAnc)++;
            continue;
        }
        if (anchorDup) {
            r.st.dupAnchor++;
            continue;
        }
        if (!dups.empty()) r.st.masked++;
        auto keep = [&](int32_t row) {
            if (row < 0) return false;
            for (const auto &m : masked)
                if (row >= m.first && row < m.second) return false;
            return true;
        };
        vector<pair<int32_t, char>> b[3];
        for (int k = 0; k < 3; k++) {
            for (const Entry &e : col(ix[k]))
                if (keep(cx.leafOf[e.genome])) b[k].push_back(make_pair(cx.leafOf[e.genome], e.base));
            sort(b[k].begin(), b[k].end());
        }
        bool same = b[0].size() == b[1].size() && b[0].size() == b[2].size();
        for (size_t j = 0; j < b[0].size() && same; j++) same = b[0][j].first == b[1][j].first && b[0][j].first == b[2][j].first;
        if (!same) {
            r.st.gap++;
            continue;
        }
        /* in the codon's orientation: the anchor's prefix, and every leaf's */
        const char a0 = oriented(upperBase(adna[(size_t)(c->p[0] - lo)])), a1 = oriented(upperBase(adna[(size_t)(c->p[1] - lo)]));
        bool agree = true;
        for (size_t j = 0; j < b[0].size() && agree; j++) {
            char x = oriented(b[0][j].second), y = oriented(b[1][j].second);
            agree = (x == 'N' || x == a0) && (y == 'N' || y == a1);
        }
        if (!agree) {
            r.st.differs++;
            continue;
        }
        r.st.passed++;
        if (b[2].size() == 1) r.st.alone++;
        Kept k;
        k.col = k.site = top[ix[2]];
        k.minus = c->minus != (topRev[ix[2]] != 0);
        k.hash = allHash[allAt(c->p[2])];
        k.tuple.assign(cx.leaves, '*');
        for (const auto &e : b[2]) k.tuple[(size_t)e.first] = oriented(e.second);
        r.kept.push_back(move(k));
    }
    r.st.testSecs += chrono::duration<double>(chrono::steady_clock::now() - t2).count();
    return r;
}

uint64_t fnv(const string &s, uint64_t h) {
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return h;
}
}

int main_4d(int argc, char **argv) {
    return runCommand("4d", [&]() {
        string sitesPath, statsPath, dupStatsPath;
        Context::Dupes dupes = Context::SPLIT;
        bool allParts = false;
        int64_t sample = -1, threads = 1;
        vector<string> exclude;
        enum { SAMPLE = 1000, SITES, STATS, EXCLUDE_GENOMES, THREADS, DUPES, DUP_STATS, PARTS };
        static const option longOpts[] = {{"sample", required_argument, nullptr, SAMPLE},
                                          {"sites", required_argument, nullptr, SITES},
                                          {"stats", required_argument, nullptr, STATS},
                                          {"excludeGenomes", required_argument, nullptr, EXCLUDE_GENOMES},
                                          {"threads", required_argument, nullptr, THREADS},
                                          {"dupes", required_argument, nullptr, DUPES},
                                          {"dupStats", required_argument, nullptr, DUP_STATS},
                                          {"parts", required_argument, nullptr, PARTS},
                                          {"help", no_argument, nullptr, 'h'},
                                          {nullptr, 0, nullptr, 0}};
        int c;
        optind = 1;
        while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
            switch (c) {
            case SAMPLE: sample = parseInt("sample", optarg); break;
            case SITES: sitesPath = optarg; break;
            case STATS: statsPath = optarg; break;
            case EXCLUDE_GENOMES: {
                stringstream ss(optarg);
                string x;
                while (getline(ss, x, ','))
                    if (!x.empty()) exclude.push_back(x);
                break;
            }
            case THREADS: threads = parseInt("threads", optarg); break;
            case DUPES:
                if (string(optarg) == "split") dupes = Context::SPLIT;
                else if (string(optarg) == "mask") dupes = Context::MASK;
                else if (string(optarg) == "drop") dupes = Context::DROP;
                else throw CliError("--dupes must be split, mask or drop");
                break;
            case DUP_STATS: dupStatsPath = optarg; break;
            case PARTS:
                if (string(optarg) == "largest") allParts = false;
                else if (string(optarg) == "all") allParts = true;
                else throw CliError("--parts must be largest or all");
                break;
            case 'h': cout << FOURD_USAGE; return 0;
            default: throw CliError("bad option");
            }
        }
        if (argc - optind < 3) throw CliError("expected an alignment, an output and at least one genome=annotation");
        if (threads < 1 || (sample != -1 && sample < 1)) throw CliError("--threads and --sample must be positive");
        Alignment aln(argv[optind]);
        const string outPath = argv[optind + 1];

        Context cx;
        cx.aln = &aln;
        cx.dupes = dupes;
        cx.allParts = allParts;
        cx.walkOpts.exclude = exclude;   /* dupes and ancestors: the defaults */
        set<const Genome *> excluded;
        for (const string &x : exclude) {
            const Genome *g = aln.genome(x);
            if (!g) throw runtime_error("no genome " + x);
            if (!g->children().empty()) throw CliError("only leaves can be excluded (" + x + " is an ancestor)");
            excluded.insert(g);
        }
        /* genomes in preorder; the leaves' rows in that order */
        vector<const Genome *> genomes, stack{aln.root()};
        while (!stack.empty()) {
            const Genome *g = stack.back();
            stack.pop_back();
            cx.index[g] = (uint32_t)genomes.size();
            genomes.push_back(g);
            for (auto ch = g->children().rbegin(); ch != g->children().rend(); ++ch) stack.push_back(*ch);
        }
        cx.genomes = genomes;
        cx.excluded.assign(genomes.size(), 0);
        cx.depth.assign(genomes.size(), 0);
        for (size_t i = 0; i < genomes.size(); i++) {
            cx.excluded[i] = excluded.count(genomes[i]) ? 1 : 0;
            cx.depth[i] = genomes[i]->depth();
        }
        vector<string> names;
        map<string, string> phastNames;
        cx.leafOf.assign(genomes.size(), -1);
        for (size_t i = 0; i < genomes.size(); i++) {
            const Genome *g = genomes[i];
            if (!g->children().empty() || excluded.count(g)) continue;
            string pn = g->name().substr(0, g->name().find('.'));
            if (phastNames.count(pn)) throw runtime_error("leaves " + phastNames[pn] + " and " + g->name() + " have the same phast name " + pn);
            phastNames[pn] = g->name();
            cx.leafOf[i] = (int32_t)names.size();
            names.push_back(pn);
        }
        cx.leaves = names.size();
        /* each genome's leaves' rows: contiguous, as rows are in preorder */
        cx.rows.assign(genomes.size(), make_pair(INT32_MAX, INT32_MIN));
        for (size_t i = genomes.size(); i-- > 0;) {
            if (cx.leafOf[i] >= 0) cx.rows[i] = make_pair(cx.leafOf[i], cx.leafOf[i] + 1);
            const Genome *p = genomes[i]->parent();
            if (p && cx.rows[i].first < cx.rows[i].second) {
                pair<int32_t, int32_t> &q = cx.rows[cx.index.at(p)];
                q.first = min(q.first, cx.rows[i].first);
                q.second = max(q.second, cx.rows[i].second);
            }
        }
        vector<uint64_t> dupBy(genomes.size(), 0);

        vector<pair<const Genome *, string>> anns;
        for (int i = optind + 2; i < argc; i++) {
            string a = argv[i];
            size_t eq = a.find('=');
            if (eq == string::npos) throw CliError("expected genome=annotation, got " + a);
            const Genome *g = aln.genome(a.substr(0, eq));
            if (!g) throw runtime_error("no genome " + a.substr(0, eq));
            if (!g->children().empty()) throw CliError(g->name() + " is an ancestor: annotations are of leaves");
            if (excluded.count(g)) throw CliError(g->name() + " is excluded");
            anns.push_back(make_pair(g, a.substr(eq + 1)));
        }

        /* the union: sites by their column, and the roles annotations give columns */
        struct Site {
            ColKey key, col;   /* the site; its third position's column (for roles) */
            bool minus;
            uint32_t first, last, nann;
            string tuple;
            uint64_t hash;
            bool live;
        };
        vector<Site> sites;
        vector<uint32_t> freeSlots;
        unordered_map<ColKey, uint32_t, ColKeyHash> siteOf;
        /* --sample N: only the 2N sites that hash lowest are held (a site that hashes higher
         * than the 2N-th can never be sampled); the N are taken from them, after conflicts */
        const size_t cap = sample > 0 ? 2 * (size_t)sample : SIZE_MAX;
        set<pair<uint64_t, uint32_t>> byHash;
        uint64_t discarded = 0;

        struct RoleOf {
            uint64_t hash;   /* the column's */
            uint8_t bits;
        };
        unordered_map<ColKey, RoleOf, ColKeyHash> roles;
        vector<pair<ColKey, uint64_t>> tested;   /* this annotation's tested site columns (cx.tried from the next) */
        vector<AnnStats> stats(anns.size());
        vector<Worker> workers((size_t)threads);
        for (size_t ai = 0; ai < anns.size(); ai++) {
            const Genome *anchor = anns[ai].first;
            AnnStats &st = stats[ai];
            vector<Codon> codons = readCodons(anns[ai].second, anchor, st);
            /* what the workers see of the union, fixed while they work on this annotation (so
             * what they do depends on no thread's timing): the cutoff and the sites held.  Roles
             * of columns that can no longer be sites go. */
            cx.cutoff = byHash.size() >= cap ? byHash.rbegin()->first : UINT64_MAX;
            cx.held.clear();
            if (!allParts) {
                for (auto &so : siteOf) cx.held.insert(so.first);
                for (auto &t : tested)
                    if (t.second <= cx.cutoff) cx.tried.insert(t);
                vector<pair<ColKey, uint64_t>>().swap(tested);
                if (cx.cutoff != UINT64_MAX) {
                    for (auto it = roles.begin(); it != roles.end();) it = it->second.hash > cx.cutoff ? roles.erase(it) : next(it);
                    for (auto it = cx.tried.begin(); it != cx.tried.end();) it = it->second > cx.cutoff ? cx.tried.erase(it) : next(it);
                }
            }
            /* units: a sequence's codons, a megabase or 1000 codons at a time */
            vector<Unit> units;
            for (size_t i = 0; i < codons.size();) {
                size_t j = i + 1;
                while (j < codons.size() && codons[j].seq == codons[i].seq && codons[j].p[0] - codons[i].p[0] < 1000000 && j - i < 1000) j++;
                units.push_back(Unit{i, j});
                i = j;
            }
            /* batches of consecutive units on one thread: its cursors and DNA handles stay warm */
            const size_t batch = max<size_t>(1, min<size_t>(16, units.size() / (32 * (size_t)threads))),
                         nb = (units.size() + batch - 1) / batch;
            size_t consumed = 0;
            const auto started = chrono::steady_clock::now();
            orderedParallel<UnitResult>(
                nb, (int)threads,
                [&](size_t b, int w) {
                    Worker &wk = workers[(size_t)w];
                    if (wk.dna) wk.dna->recycle();
                    UnitResult all;
                    for (size_t k = b * batch; k < min(units.size(), (b + 1) * batch); k++) {
                        UnitResult r = runUnit(cx, wk, anchor, codons, units[k]);
                        all.st.add(r.st);
                        for (auto &d : r.dupBy) all.dupBy[d.first] += d.second;
                        for (Kept &x : r.kept) all.kept.push_back(move(x));
                        all.roles.insert(all.roles.end(), r.roles.begin(), r.roles.end());
                        all.tested.insert(all.tested.end(), r.tested.begin(), r.tested.end());
                    }
                    return all;
                },
                [&](size_t, UnitResult &r) {
                    st.add(r.st);
                    for (auto &d : r.dupBy) dupBy[d.first] += d.second;
                    for (auto &t : r.tested)
                        if (byHash.size() < cap || t.second <= byHash.rbegin()->first) tested.push_back(t);
                    for (const Role &x : r.roles) {
                        RoleOf &ro = roles[x.col];
                        ro.hash = x.hash;
                        ro.bits |= x.bit;
                    }
                    for (Kept &k : r.kept) {
                        auto it = siteOf.find(k.site);
                        if (it != siteOf.end()) {
                            if (sites[it->second].last != ai) {
                                sites[it->second].last = (uint32_t)ai;
                                sites[it->second].nann++;
                            }
                            continue;
                        }
                        if (k.tuple.empty()) continue;   /* held as the annotation started, since gone */
                        const uint64_t h = k.hash;
                        if (byHash.size() >= cap && h >= byHash.rbegin()->first) {
                            discarded++;
                            continue;
                        }
                        uint32_t slot;
                        if (!freeSlots.empty()) {
                            slot = freeSlots.back();
                            freeSlots.pop_back();
                        } else {
                            slot = (uint32_t)sites.size();
                            sites.push_back(Site());
                        }
                        sites[slot] = Site{k.site, k.col, k.minus, (uint32_t)ai, (uint32_t)ai, 1, move(k.tuple), h, true};
                        siteOf[k.site] = slot;
                        byHash.insert(make_pair(h, slot));
                        if (byHash.size() > cap) {   /* the highest goes */
                            const uint32_t e = byHash.rbegin()->second;
                            byHash.erase(prev(byHash.end()));
                            siteOf.erase(sites[e].key);
                            sites[e].live = false;
                            string().swap(sites[e].tuple);
                            freeSlots.push_back(e);
                            discarded++;
                        }
                    }
                    consumed++;
                    if (consumed * 100 / nb != (consumed - 1) * 100 / nb)
                        cerr << "hal 4d: " << anchor->name() << ": " << consumed * 100 / nb << "% ("
                             << (int)chrono::duration<double>(chrono::steady_clock::now() - started).count() << " s; widest column "
                             << st.widest << " bases)" << endl;
                });
            cerr << "hal 4d: " << anchor->name() << ": " << st.transcripts << " transcripts, " << st.codons << " codons, " << st.candidates
                 << " fourfold degenerate there (" << st.held << " at sites already held, " << st.tried << " at columns already tested, " << st.above
                 << " above the sample), " << st.passed << " kept; " << siteOf.size() << " sites held ("
                 << (int)chrono::duration<double>(chrono::steady_clock::now() - started).count() << " s; thread time: lift " << (int)st.liftSecs
                 << " s, columns " << (int)st.sweepSecs << " s, tests " << (int)st.testSecs << " s)" << endl;
        }

        /* drop columns given two roles; then the sample: the lowest hashes of the columns */
        vector<uint8_t> status(sites.size(), 0);   /* 0 kept, 1 conflict, 2 not sampled */
        uint64_t conflicts = 0;
        vector<pair<uint64_t, uint32_t>> order;
        for (size_t i = 0; i < sites.size(); i++) {
            const Site &s = sites[i];
            if (!s.live) continue;
            auto ro = roles.find(s.col);
            if (ro == roles.end() || ro->second.bits != (uint8_t)(1 << (4 + (s.minus ? 1 : 0)))) {
                status[i] = 1;
                conflicts++;
                continue;
            }
            order.push_back(make_pair(s.hash, (uint32_t)i));
        }
        sort(order.begin(), order.end());
        if (sample > 0 && (size_t)sample < order.size()) {
            for (size_t i = (size_t)sample; i < order.size(); i++) status[order[i].second] = 2;
            order.resize((size_t)sample);
        }
        /* the SS: kept sites by column, identical tuples counted together */
        vector<uint32_t> keep;
        for (auto &o : order) keep.push_back(o.second);
        sort(keep.begin(), keep.end(), [&](uint32_t a, uint32_t b) { return sites[a].key < sites[b].key; });
        vector<pair<const string *, uint64_t>> tuples;
        unordered_map<string, size_t> tupleOf;
        for (uint32_t i : keep) {
            auto it = tupleOf.find(sites[i].tuple);
            if (it == tupleOf.end()) {
                tupleOf[sites[i].tuple] = tuples.size();
                tuples.push_back(make_pair(&sites[i].tuple, 1));
            } else {
                tuples[it->second].second++;
            }
        }
        ofstream out(outPath);
        if (!out) throw runtime_error("cannot write " + outPath);
        out << "NSEQS = " << names.size() << "\nLENGTH = " << keep.size() << "\nTUPLE_SIZE = 1\nNTUPLES = " << tuples.size() << "\nNAMES = ";
        for (size_t i = 0; i < names.size(); i++) out << (i ? "," : "") << names[i];
        out << "\nALPHABET = ACGT\nNCATS = -1\n\n";
        for (size_t i = 0; i < tuples.size(); i++) out << i << '\t' << *tuples[i].first << '\t' << tuples[i].second << '\n';
        out.close();
        if (!out) throw runtime_error("cannot write " + outPath);

        if (!sitesPath.empty()) {
            ofstream s(sitesPath);
            if (!s) throw runtime_error("cannot write " + sitesPath);
            s << "#top\tpos\tstrand\tannotations\tfirst\tstatus\n";
            static const char *what[] = {"kept", "conflict", "unsampled"};
            vector<uint32_t> all;
            for (size_t i = 0; i < sites.size(); i++)
                if (sites[i].live) all.push_back((uint32_t)i);
            sort(all.begin(), all.end(), [&](uint32_t a, uint32_t b) { return sites[a].key < sites[b].key; });
            for (uint32_t i : all) {
                const Site &x = sites[i];
                const Genome *g = genomes[x.key.genome];
                s << g->name() << '.' << g->sequenceName(x.key.seq) << '\t' << x.key.pos << '\t' << (x.minus ? '-' : '+') << '\t' << x.nann << '\t'
                  << anns[x.first].first->name() << '\t' << what[status[i]] << '\n';
            }
            if (!s) throw runtime_error("cannot write " + sitesPath);
        }
        if (!statsPath.empty()) {
            ofstream s(statsPath);
            if (!s) throw runtime_error("cannot write " + statsPath);
            s << "genome\tannotation\tcdsLines\ttranscripts\tunknownSequence\tphaseDisagrees\tcodons\tframeConflict\tanchorN\tnotFourfold\t"
                 "fourfold\talreadyHeld\talreadyTested\taboveSample\tdupLeaf\tdupAncestor\tdupAnnotated\tmaskedDup\tspliced\tparts\tpartsAlone\tsplitColumns\tmissingBase\tprefixDiffers\tkept\tkeptAlone\t"
                 "liftSecs\tsweepSecs\ttestSecs\n";
            for (size_t i = 0; i < anns.size(); i++) {
                const AnnStats &t = stats[i];
                s << anns[i].first->name() << '\t' << anns[i].second << '\t' << t.cdsLines << '\t' << t.transcripts << '\t' << t.unknownSeq << '\t'
                  << t.phaseBad << '\t' << t.codons << '\t' << t.frameConflict << '\t' << t.anchorN << '\t' << t.notPrefix << '\t' << t.candidates
                  << '\t' << t.held << '\t' << t.tried << '\t' << t.above << '\t' << t.dupLeaf << '\t' << t.dupAnc << '\t' << t.dupAnchor << '\t' << t.masked << '\t' << t.spliced << '\t' << t.parts << '\t' << t.partsAlone << '\t' << t.columns << '\t' << t.gap << '\t' << t.differs << '\t' << t.passed << '\t' << t.alone << '\t'
                  << t.liftSecs << '\t' << t.sweepSecs << '\t' << t.testSecs << '\n';
            }
            if (!s) throw runtime_error("cannot write " + statsPath);
        }
        if (!dupStatsPath.empty()) {
            ofstream s(dupStatsPath);
            if (!s) throw runtime_error("cannot write " + dupStatsPath);
            s << "genome\tleaf\tleaves\ttestedCodonsWithTwoCopies\n";
            vector<size_t> idx;
            for (size_t i = 0; i < genomes.size(); i++)
                if (dupBy[i]) idx.push_back(i);
            sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return dupBy[a] > dupBy[b] || (dupBy[a] == dupBy[b] && a < b); });
            for (size_t i : idx)
                s << genomes[i]->name() << '\t' << (genomes[i]->children().empty() ? "leaf" : "ancestor") << '\t'
                  << max(0, cx.rows[i].second - cx.rows[i].first) << '\t' << dupBy[i] << '\n';
            if (!s) throw runtime_error("cannot write " + dupStatsPath);
        }
        cerr << "hal 4d: " << siteOf.size() << " sites held" << (sample > 0 ? " (the lowest hashes; " + to_string(discarded) + " passed over)" : "")
             << ", " << conflicts << " of them given two roles (dropped), " << keep.size() << " written (" << tuples.size() << " distinct)" << endl;
        if (sample > 0 && keep.size() < (size_t)sample && discarded > 0)
            cerr << "hal 4d: warning: " << keep.size() << " sites, not " << sample << ": of the " << cap << " held, more than half were given two roles"
                 << endl;
        return 0;
    });
}
}
