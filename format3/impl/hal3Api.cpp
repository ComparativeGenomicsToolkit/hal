/*
 * HAL format 3 API: the public classes (inc/hal3/), over the internal ones.
 */
#include "hal3/hal3.h"

#include "hal3Alignment.h"
#include "hal3Liftover.h"
#include "hal3Walk.h"

#include <algorithm>
#include <cstring>
#include <mutex>

using namespace std;

namespace hal3 {
namespace api {

/* ---------------------------------------------------------------- the alignment */

struct Alignment::Impl {
    unique_ptr<hal3::Alignment> aln;
    vector<unique_ptr<Genome>> owned;
    vector<const Genome *> genomes;
    map<string, const Genome *> byName;
    map<const hal3::Genome *, const Genome *> of;   /* internal -> public */
    const Genome *pub(const hal3::Genome *g) const {
        auto i = of.find(g);
        return i == of.end() ? nullptr : i->second;
    }
};

struct Genome::Impl {
    const hal3::Genome *g = nullptr;
    const GenomeEntry *entry = nullptr;
    const Genome *parent = nullptr;
    vector<const Genome *> children;
    mutable once_flag named;
    mutable vector<string> names;
    const vector<string> &sequenceNames() const {
        call_once(named, [this]() {   /* nothing thrown out of call_once (see hal3Alignment.cpp) */
            try {
                uint32_t n = g->numSequences();
                names.reserve(n);
                for (uint32_t i = 0; i < n; i++) names.push_back(g->sequenceName(i));
            } catch (...) {
                names.clear();
            }
        });
        if (names.size() != (size_t)entry->numSequences) {
            g->numSequences();   /* rethrows why the file could not be read, if it could not */
            throw Error(ErrorKind::Format, "genome " + entry->name + ": the manifest says " + to_string(entry->numSequences) +
                                               " sequences, its genome file " + to_string(names.size()));
        }
        return names;
    }
    void checkSeq(uint32_t seq) const {
        if (seq >= entry->numSequences)
            throw Error(ErrorKind::Argument, "genome " + entry->name + " has no sequence " + to_string(seq));
    }
};

Alignment::Alignment() : _impl(new Impl()) {}
Alignment::~Alignment() {}

shared_ptr<const Alignment> Alignment::open(const string &path) {
    shared_ptr<Alignment> a(new Alignment());
    Impl &m = *a->_impl;
    m.aln.reset(new hal3::Alignment(path));
    for (const GenomeEntry &e : m.aln->manifest().genomes) {
        const hal3::Genome *g = m.aln->genome(e.name);
        unique_ptr<Genome> p(new Genome());
        p->_impl->g = g;
        p->_impl->entry = &e;
        m.of[g] = p.get();
        m.byName[e.name] = p.get();
        m.genomes.push_back(p.get());
        m.owned.push_back(move(p));
    }
    for (auto &p : m.owned) {
        Genome::Impl &gi = *p->_impl;
        gi.parent = m.pub(gi.g->parent());
        for (const hal3::Genome *c : gi.g->children()) gi.children.push_back(m.pub(c));
    }
    return a;
}

const string &Alignment::path() const { return _impl->aln->path(); }
bool Alignment::isPack() const { return _impl->aln->store().isPack(); }
const Genome *Alignment::root() const { return _impl->pub(_impl->aln->root()); }
const Genome *Alignment::genome(const string &name) const {
    auto i = _impl->byName.find(name);
    return i == _impl->byName.end() ? nullptr : i->second;
}
const vector<const Genome *> &Alignment::genomes() const { return _impl->genomes; }
const string &Alignment::newick() const { return _impl->aln->manifest().newick; }
const map<string, string> &Alignment::meta() const { return _impl->aln->manifest().meta; }

/* ---------------------------------------------------------------- genomes */

Genome::Genome() : _impl(new Impl()) {}
Genome::~Genome() {}

const string &Genome::name() const { return _impl->entry->name; }
const Genome *Genome::parent() const { return _impl->parent; }
const vector<const Genome *> &Genome::children() const { return _impl->children; }
int Genome::depth() const { return _impl->g->depth(); }
double Genome::branchLength() const { return _impl->entry->branchLength; }
int64_t Genome::length() const { return (int64_t)_impl->entry->length; }
int64_t Genome::numSequences() const { return (int64_t)_impl->entry->numSequences; }
bool Genome::hasDna() const { return !_impl->entry->dnaPath.empty(); }

vector<string> Genome::memberRoles() const {
    vector<string> r;
    for (const MemberRef &m : _impl->entry->members) r.push_back(m.role);
    return r;
}

const string &Genome::sequenceName(uint32_t seq) const {
    _impl->checkSeq(seq);
    return _impl->sequenceNames()[seq];
}

int64_t Genome::sequenceLength(uint32_t seq) const {
    _impl->checkSeq(seq);
    return _impl->g->sequenceLength(seq);
}

int64_t Genome::sequenceStart(uint32_t seq) const {
    _impl->checkSeq(seq);
    return _impl->g->sequenceStart(seq);
}

int64_t Genome::findSequence(const string &name) const { return _impl->g->findSequence(name); }

uint32_t Genome::sequenceAt(int64_t genomePos) const {
    if (genomePos < 0 || genomePos >= length())
        throw Error(ErrorKind::Argument, "genome " + name() + " has no position " + to_string(genomePos));
    uint32_t lo = 0, hi = (uint32_t)numSequences();   /* the last sequence starting at or before genomePos */
    while (hi - lo > 1) {
        uint32_t mid = lo + (hi - lo) / 2;
        if (_impl->g->sequenceStart(mid) <= genomePos) lo = mid;
        else hi = mid;
    }
    return lo;
}

map<string, string> Genome::meta() const {
    map<string, string> m;
    lmap_file *f = _impl->g->hgn();
    const string prefix = "hal.meta.";
    for (uint32_t i = 0; i < lmap_meta_count(f); i++) {
        const char *k, *v;
        if (lmap_meta_at(f, i, &k, &v) == 0 && !strncmp(k, prefix.c_str(), prefix.size())) m[k + prefix.size()] = v;
    }
    return m;
}

const Genome *lowestCommonAncestor(const Genome *a, const Genome *b) {
    while (a->depth() > b->depth()) a = a->parent();
    while (b->depth() > a->depth()) b = b->parent();
    while (a != b) {
        a = a->parent();
        b = b->parent();
    }
    return a;
}

/* ---------------------------------------------------------------- liftover */

namespace {

const hal3::Genome *internal(const Genome *g) {
    if (!g) throw Error(ErrorKind::Argument, "no genome (nullptr)");
    return g->impl().g;
}

void toInternal(const vector<LiftPiece> &in, vector<hal3::LiftPiece> &out) {
    out.clear();
    out.reserve(in.size());
    for (const LiftPiece &p : in) out.push_back(hal3::LiftPiece{p.srcStart, p.tgtSeq, p.tgtStart, p.length, p.reversed});
}

void toPublic(const vector<hal3::LiftPiece> &in, vector<LiftPiece> &out) {
    out.clear();
    out.reserve(in.size());
    for (const hal3::LiftPiece &p : in) out.push_back(LiftPiece{p.srcStart, p.tgtSeq, p.tgtStart, p.length, p.reversed});
}
}

struct Lifter::Impl {
    shared_ptr<const Alignment> aln;
    CursorCache cursors;
    vector<hal3::LiftPiece> a, b, c;
    explicit Impl(shared_ptr<const Alignment> al) : aln(move(al)) {}
};

Lifter::Lifter(shared_ptr<const Alignment> alignment) : _impl(new Impl(move(alignment))) {}
Lifter::~Lifter() {}

void Lifter::lift(const Genome *src, uint32_t srcSeq, int64_t start, int64_t end, const Genome *tgt, bool dupes,
                  vector<LiftPiece> &out) {
    const hal3::Genome *s = internal(src), *t = internal(tgt);
    if (srcSeq >= (uint64_t)src->numSequences() || start < 0 || end < start || end > src->sequenceLength(srcSeq))
        throw Error(ErrorKind::Argument, "no interval [" + to_string(start) + ", " + to_string(end) + ") of sequence " +
                                             to_string(srcSeq) + " of " + src->name());
    hal3::liftInterval(s, srcSeq, start, end, t, dupes, _impl->a, &_impl->cursors);
    toPublic(_impl->a, out);
}

void Lifter::up(const Genome *g, const vector<LiftPiece> &in, vector<LiftPiece> &out, vector<LiftPiece> *unmapped,
                bool canonicalOnly) {
    toInternal(in, _impl->a);
    hal3::liftUpBranch(internal(g), _impl->a, _impl->b, unmapped ? &_impl->c : nullptr, _impl->cursors, canonicalOnly);
    toPublic(_impl->b, out);
    if (unmapped) toPublic(_impl->c, *unmapped);
}

void Lifter::down(const Genome *child, bool dupes, const vector<LiftPiece> &in, vector<LiftPiece> &out) {
    toInternal(in, _impl->a);
    hal3::liftDownBranch(internal(child), dupes, _impl->a, _impl->b, _impl->cursors);
    toPublic(_impl->b, out);
}

/* ---------------------------------------------------------------- the reference walk */

struct ReferenceWalk::Impl {
    shared_ptr<const Alignment> aln;
    unique_ptr<TreeWalk> walk;
    set<const Genome *> printed;
    PieceMap lifted;
};

ReferenceWalk::ReferenceWalk(shared_ptr<const Alignment> alignment, const Genome *ref, const WalkOptions &o)
    : _impl(new Impl()) {
    _impl->aln = move(alignment);
    hal3::WalkOptions w;
    w.ancestors = o.ancestors;
    w.dupes = o.dupes;
    w.onlyOrthologs = o.onlyOrthologs;
    w.targets = o.targets;
    w.rootGenome = o.rootGenome;
    w.exclude = o.exclude;
    w.refSingleCopy = o.refSingleCopy;
    const hal3::Alignment &aln = *_impl->aln->impl().aln;
    try {
        _impl->walk.reset(new TreeWalk(aln, internal(ref), w));
    } catch (Error &) {
        throw;
    } catch (runtime_error &e) {
        throw Error(ErrorKind::Argument, e.what());
    }
    for (const hal3::Genome *g : _impl->walk->printed()) _impl->printed.insert(_impl->aln->impl().pub(g));
}

ReferenceWalk::~ReferenceWalk() {}

void ReferenceWalk::lift(uint32_t refSeq, const vector<pair<int64_t, int64_t>> &intervals,
                         map<const Genome *, vector<LiftPiece>> &pieces) {
    const Genome *ref = _impl->aln->impl().pub(_impl->walk->ref());
    if (refSeq >= (uint64_t)ref->numSequences())
        throw Error(ErrorKind::Argument, "genome " + ref->name() + " has no sequence " + to_string(refSeq));
    int64_t prev = 0, len = ref->sequenceLength(refSeq);
    for (const auto &iv : intervals)
        if (iv.first < prev || iv.second < iv.first || iv.second > len)
            throw Error(ErrorKind::Argument, "reference intervals must be sorted, disjoint and inside the sequence");
        else prev = iv.second;
    _impl->walk->lift(refSeq, intervals, _impl->lifted);
    pieces.clear();
    for (auto &gp : _impl->lifted) toPublic(gp.second, pieces[_impl->aln->impl().pub(gp.first)]);
}

const set<const Genome *> &ReferenceWalk::printed() const { return _impl->printed; }
}
}
