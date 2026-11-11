/*
 * HAL format 3: universal columns as trees of bases.  See hal3Columns.h.
 */
#include "hal3Columns.h"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

using namespace std;

namespace hal3 {

struct UniversalColumns::Worker {
    TreeWalk walk;
    DnaReader dna;
    vector<Segment> segs;
    vector<int32_t> local;   /* per base of the chunk's pieces: its entry, counted from its column's first */
    vector<int64_t> fill;    /* per column: entries placed so far */
    string bases;
    Worker(const Alignment &aln, const Genome *backbone, const WalkOptions &o) : walk(aln, backbone, o), dna(aln) {}
};

UniversalColumns::UniversalColumns(const Alignment &aln, const vector<RefInterval> &ranges, const ColumnOptions &opts,
                                   int workers)
    : _aln(aln), _opts(opts) {
    if (workers < 1) throw runtime_error("UniversalColumns: workers must be positive");
    if (opts.chunk < 1) throw runtime_error("UniversalColumns: the chunk must be positive");
    const Genome *backbone = opts.rootGenome.empty() ? aln.root() : aln.genome(opts.rootGenome);
    if (!backbone) throw runtime_error("no genome " + opts.rootGenome);
    set<const Genome *> excluded;
    for (const string &x : opts.exclude) {
        const Genome *g = aln.genome(x);
        if (!g) throw runtime_error("no genome " + x);
        if (!g->children().empty()) throw runtime_error("only leaves can be excluded (" + x + " is an ancestor)");
        if (g == backbone) throw runtime_error(x + " is the backbone: it cannot be excluded");
        excluded.insert(g);
    }
    /* the genomes, in preorder from the backbone */
    vector<const Genome *> stack{backbone};
    while (!stack.empty()) {
        const Genome *g = stack.back();
        stack.pop_back();
        if (excluded.count(g)) continue;
        _indexOf[g] = (uint32_t)_genomes.size();
        _genomes.push_back(g);
        _parentOf.push_back(g == backbone ? -1 : (int32_t)_indexOf.at(g->parent()));
        for (auto c = g->children().rbegin(); c != g->children().rend(); ++c) stack.push_back(*c);
    }
    _walkOpts.universal = true;
    _walkOpts.dupes = true;
    _walkOpts.ancestors = true;
    _walkOpts.allAncestorRows = true;
    _walkOpts.rootGenome = opts.rootGenome;
    _walkOpts.exclude = opts.exclude;
    { TreeWalk check(aln, backbone, _walkOpts); }   /* bad options fail before the index is built */
    const bool whole = ranges.empty();
    _index.reset(new UniversalIndex(aln, backbone, opts.leafInsertions, true, excluded, whole, workers));
    vector<TopRange> top;
    if (whole) {
        for (uint32_t i = 0; i < backbone->numSequences(); i++) top.push_back(TopRange{0, i, 0, backbone->sequenceLength(i)});
    } else {
        for (const RefInterval &iv : ranges) top.push_back(TopRange{0, iv.seq, iv.start, iv.end});
    }
    _plan.reset(new UniversalPlanner(*_index, top, whole, opts.chunk));
    _plan->run();
    _axes = axisRanges(*_index, *_plan, opts.chunk);
    _workers.resize((size_t)workers);
}

UniversalColumns::~UniversalColumns() {}

void UniversalColumns::chunk(size_t k, int w, ColumnChunk &out) {
    if (k >= chunks()) throw runtime_error("UniversalColumns: no chunk " + to_string(k));
    if (w < 0 || w >= (int)_workers.size()) throw runtime_error("UniversalColumns: no worker " + to_string(w));
    unique_ptr<Worker> &wp = _workers[(size_t)w];
    if (!wp) wp.reset(new Worker(_aln, backbone(), _walkOpts));
    Worker &wk = *wp;
    wk.dna.recycle();

    /* the chunk's tops: placed intervals (own) and the top-level range's bases (quiet),
     * then every base below them, each piece with the parent piece it came from */
    _plan->segments(k, wk.segs);
    const UniversalPlanner::Chunk &c = _plan->chunk(k);
    PieceMap own, quiet, pieces;
    map<const Genome *, vector<int32_t>> from;
    for (const Segment &sg : wk.segs)
        (sg.depth ? own : quiet)[_index->genome(sg.g)].push_back(LiftPiece{sg.u, sg.seq, sg.pos, sg.len, sg.rev});
    wk.walk.downParented(own, quiet, pieces, from);

    const AxisRange &a = _axes[c.range];
    out.axisName = &a.name;
    out.axisStart = a.offset + (c.u0 - a.u0);
    out.cols = c.cols;

    /* the genomes with pieces, in preorder; each piece's first base in wk.local */
    vector<pair<uint32_t, const Genome *>> order;
    for (auto &gp : pieces)
        if (!gp.second.empty() && _indexOf.count(gp.first)) order.push_back(make_pair(_indexOf.at(gp.first), gp.first));
    sort(order.begin(), order.end());
    map<const Genome *, vector<int64_t>> at;   /* per piece: its first base in wk.local */
    int64_t nbases = 0;
    out.colStart.assign((size_t)c.cols + 1, 0);
    for (auto &og : order) {
        vector<int64_t> &v = at[og.second];
        for (const LiftPiece &p : pieces.at(og.second)) {
            if (p.srcStart < 0 || p.srcStart + p.length > c.cols)
                throw runtime_error("internal error: a piece of " + og.second->name() + " is outside its chunk");
            v.push_back(nbases);
            nbases += p.length;
            for (int64_t u = p.srcStart; u < p.srcStart + p.length; u++) out.colStart[(size_t)u + 1]++;
        }
    }
    for (int64_t u = 0; u < c.cols; u++) out.colStart[(size_t)u + 1] += out.colStart[(size_t)u];
    out.genome.assign((size_t)nbases, 0);
    out.parent.assign((size_t)nbases, -1);
    out.piece.assign((size_t)nbases, 0);
    out.base.assign((size_t)nbases, 'N');
    out.pieces.clear();
    wk.local.assign((size_t)nbases, -1);
    wk.fill.assign((size_t)c.cols, 0);

    for (auto &og : order) {
        const uint32_t gi = og.first;
        const Genome *g = og.second;
        const vector<LiftPiece> &ps = pieces.at(g);
        const vector<int32_t> &fs = from.at(g);
        const vector<int64_t> &mine = at.at(g);
        const Genome *pg = g->parent();
        const vector<LiftPiece> *pps = nullptr;
        const vector<int64_t> *pat = nullptr;
        if (g != backbone()) {
            auto pp = pieces.find(pg);
            if (pp != pieces.end()) {
                pps = &pp->second;
                pat = &at.at(pg);
            }
        }
        const GenomeEntry *e = _aln.manifest().find(g->name());
        const bool dna = e && !e->dnaPath.empty() && (_opts.ancestorBases || g->children().empty());
        for (size_t i = 0; i < ps.size(); i++) {
            const LiftPiece &p = ps[i];
            const uint32_t id = (uint32_t)out.pieces.size();
            out.pieces.push_back(ColumnPiece{gi, p.tgtSeq, p.tgtStart, p.length, p.srcStart, p.reversed});
            if (dna) {
                wk.dna.fetch(g, p.tgtSeq, p.tgtStart, p.tgtStart + p.length, wk.bases);
                if (p.reversed) reverseComplement(wk.bases);
            }
            const int32_t f = fs[i];
            if (f >= 0 && (!pps || (size_t)f >= pps->size()))
                throw runtime_error("internal error: a piece of " + g->name() + " has no parent piece");
            for (int64_t d = 0; d < p.length; d++) {
                const int64_t u = p.srcStart + d;
                const int64_t slot = wk.fill[(size_t)u]++;
                const int64_t ei = out.colStart[(size_t)u] + slot;
                wk.local[(size_t)(mine[i] + d)] = (int32_t)slot;
                int32_t parent = -1;
                if (f >= 0) {
                    const LiftPiece &q = (*pps)[(size_t)f];
                    const int64_t dq = u - q.srcStart;
                    if (dq < 0 || dq >= q.length)
                        throw runtime_error("internal error: a base of " + g->name() + " is not in its parent piece's column");
                    parent = wk.local[(size_t)((*pat)[(size_t)f] + dq)];
                }
                out.genome[(size_t)ei] = gi;
                out.parent[(size_t)ei] = parent;
                out.piece[(size_t)ei] = id;
                if (dna) out.base[(size_t)ei] = wk.bases[(size_t)d];
            }
        }
    }
}
}
