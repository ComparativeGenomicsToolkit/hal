/*
 * HAL format 3: the tree walk.  See hal3Walk.h.
 */
#include "hal3Walk.h"

#include <algorithm>
#include <stdexcept>
#include <tuple>

using namespace std;

namespace hal3 {

TreeWalk::TreeWalk(const Alignment &alignment, const Genome *ref, const WalkOptions &opts)
    : _aln(alignment), _ref(ref), _opts(opts), _cursors((size_t)1 << 20, 1024) {
    /* printed genomes: the targets (or the subtree of rootGenome) and the reference;
     * scope: the spanning tree of those, whose top is where the walk turns */
    auto subtree = [](const Genome *r, vector<const Genome *> &out) {
        vector<const Genome *> stack{r};
        while (!stack.empty()) {
            const Genome *g = stack.back();
            stack.pop_back();
            out.push_back(g);
            for (const Genome *c : g->children()) stack.push_back(c);
        }
    };
    vector<const Genome *> targets;
    if (!_opts.rootGenome.empty()) {
        if (!_opts.targets.empty()) throw runtime_error("--rootGenome and --targetGenomes cannot be combined");
        const Genome *r = _aln.genome(_opts.rootGenome);
        if (!r) throw runtime_error("no genome " + _opts.rootGenome);
        subtree(r, targets);
    }
    for (const string &t : _opts.targets) {
        const Genome *g = _aln.genome(t);
        if (!g) throw runtime_error("no genome " + t);
        targets.push_back(g);
    }
    if (targets.empty()) subtree(_aln.root(), targets);
    targets.push_back(_ref);
    _top = targets[0];
    for (const Genome *t : targets) _top = lowestCommonAncestor(_top, t);
    for (const Genome *t : targets) {
        for (const Genome *g = t; g != _top; g = g->parent()) _scope.insert(g);
        _printed.insert(t);
    }
    _scope.insert(_top);
    for (const Genome *g = _ref;; g = g->parent()) {
        _refPath.push_back(g);
        if (g == _top) break;
    }
    if (!_opts.ancestors) {
        if (!_ref->children().empty()) throw runtime_error("--noAncestors needs a leaf reference genome");
        for (auto it = _printed.begin(); it != _printed.end();)
            it = (*it)->children().empty() ? next(it) : _printed.erase(it);
    }
    for (const string &x : _opts.exclude) {
        const Genome *g = _aln.genome(x);
        if (!g) throw runtime_error("no genome " + x);
        if (g == _ref) throw runtime_error("the reference genome " + x + " cannot be excluded");
        _printed.erase(g);
    }
    if (_opts.refSingleCopy) _printed.erase(_ref);
}

vector<RefInterval> chunkIntervals(const vector<RefInterval> &intervals, int64_t size) {
    vector<RefInterval> out;
    for (const RefInterval &iv : intervals)
        for (int64_t s = iv.start; s < iv.end; s += size) out.push_back(RefInterval{iv.seq, s, min(iv.end, s + size)});
    return out;
}

void coalescePieces(vector<LiftPiece> &v) {
    if (v.size() < 2) return;
    auto diag = [](const LiftPiece &p) { return p.reversed ? p.srcStart + p.tgtStart + p.length : p.tgtStart - p.srcStart; };
    sort(v.begin(), v.end(), [&](const LiftPiece &a, const LiftPiece &b) {
        return make_tuple(a.tgtSeq, a.reversed, diag(a), a.srcStart) < make_tuple(b.tgtSeq, b.reversed, diag(b), b.srcStart);
    });
    vector<LiftPiece> o;
    for (const LiftPiece &p : v) {
        if (!o.empty()) {
            LiftPiece &q = o.back();
            if (q.tgtSeq == p.tgtSeq && q.reversed == p.reversed && diag(q) == diag(p) && q.srcStart + q.length == p.srcStart) {
                if (q.reversed) q.tgtStart = p.tgtStart;
                q.length += p.length;
                continue;
            }
        }
        o.push_back(p);
    }
    v.swap(o);
}

/* At v: its own pieces (arrived from the parent, plus any tops at v), recorded if v is
 * printed, then carried into each child in scope. */
void TreeWalk::descend(const Genome *v, vector<LiftPiece> &in, map<const Genome *, vector<LiftPiece>> &tops,
                       map<const Genome *, vector<LiftPiece>> &pieces) {
    auto t = tops.find(v);
    if (t != tops.end()) {
        in.insert(in.end(), t->second.begin(), t->second.end());
        tops.erase(t);
    }
    coalescePieces(in);
    if (_printed.count(v)) {
        vector<LiftPiece> &r = pieces[v];
        r.insert(r.end(), in.begin(), in.end());
    }
    /* a child with no pieces may still lead to tops (a walk up stops short of the scope
     * top when every piece tops out lower down) */
    for (const Genome *c : v->children()) {
        if (!_scope.count(c)) continue;
        vector<LiftPiece> down;
        liftDownBranch(c, _opts.dupes, in, down, _cursors);
        if (!down.empty() || _below.count(c)) descend(c, down, tops, pieces);
    }
}

void TreeWalk::lift(uint32_t refSeq, const vector<pair<int64_t, int64_t>> &intervals,
                    map<const Genome *, vector<LiftPiece>> &pieces) {
    pieces.clear();
    _below.clear();
    vector<LiftPiece> cur, next, unmapped;
    for (const auto &iv : intervals)
        if (iv.second > iv.first) cur.push_back(LiftPiece{iv.first, refSeq, iv.first, iv.second - iv.first, false});
    if (cur.empty()) return;
    map<const Genome *, vector<LiftPiece>> tops;
    if (!_opts.onlyOrthologs) {
        /* up to the scope top; every top descends into its whole subtree */
        const Genome *a = _ref;
        for (;;) {
            if (a == _top) {
                vector<LiftPiece> &t = tops[a];
                t.insert(t.end(), cur.begin(), cur.end());
                break;
            }
            liftUpBranch(a, cur, next, &unmapped, _cursors, !_opts.dupes);
            if (!unmapped.empty()) {
                vector<LiftPiece> &t = tops[a];
                t.insert(t.end(), unmapped.begin(), unmapped.end());
            }
            cur.swap(next);
            if (cur.empty()) break;
            a = a->parent();
        }
        for (auto &t : tops)
            for (const Genome *g = t.first; g && _below.insert(g).second; g = g->parent()) {
            }
        vector<LiftPiece> none;
        descend(_top, none, tops, pieces);
    } else {
        /* liftover to the common ancestor with each genome: at each ancestor on the
         * reference's lineage, descend only into the children off the lineage; and
         * into the reference's own subtree */
        const Genome *below = nullptr;
        for (const Genome *a = _ref;; a = a->parent()) {
            if (a != _ref && _printed.count(a)) {
                vector<LiftPiece> &r = pieces[a];
                r.insert(r.end(), cur.begin(), cur.end());
            }
            for (const Genome *c : a->children()) {
                if (c == below || !_scope.count(c)) continue;
                vector<LiftPiece> down;
                liftDownBranch(c, _opts.dupes, cur, down, _cursors);
                if (!down.empty()) descend(c, down, tops, pieces);
            }
            if (a == _top) break;
            liftUpBranch(a, cur, next, nullptr, _cursors, !_opts.dupes);
            cur.swap(next);
            coalescePieces(cur);
            if (cur.empty()) break;
            below = a;
        }
    }
    /* the reference's own identity pieces */
    auto r = pieces.find(_ref);
    if (r != pieces.end()) {
        vector<LiftPiece> &v = r->second;
        v.erase(remove_if(v.begin(), v.end(),
                          [&](const LiftPiece &p) { return !p.reversed && p.tgtSeq == refSeq && p.tgtStart == p.srcStart; }),
                v.end());
    }
}
}
