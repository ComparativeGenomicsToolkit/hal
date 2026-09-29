/*
 * HAL format 3: conversion from a HAL alignment.  See hal3Convert.h.
 */
#include "hal3Convert.h"
#include "hal.h"
#include "hal3/hal3Writer.h"
#include "hal3Manifest.h"
#include "hal3Store.h"
#include "liftmap.h"
#include "twobit64.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace hal;

namespace hal3 {

namespace {
struct SeqInfo {
    string name;
    hal_index_t start;     /* in genome coordinates */
    hal_size_t length;
};

vector<SeqInfo> sequencesOf(const Genome *genome) {
    vector<SeqInfo> seqs;
    for (SequenceIteratorPtr si = genome->getSequenceIterator(); !si->atEnd(); si->toNext()) {
        const Sequence *s = si->getSequence();
        seqs.push_back({s->getName(), s->getStartPosition(), s->getSequenceLength()});
    }
    return seqs;
}

/* The sequence containing genome position pos (sequences sorted by start). */
size_t sequenceAt(const vector<SeqInfo> &seqs, hal_index_t pos) {
    auto it = upper_bound(seqs.begin(), seqs.end(), pos,
                          [](hal_index_t p, const SeqInfo &s) { return p < s.start; });
    if (it == seqs.begin()) throw runtime_error("position " + to_string(pos) + " before the first sequence");
    size_t i = (size_t)(it - seqs.begin()) - 1;
    if (pos >= seqs[i].start + (hal_index_t)seqs[i].length)
        throw runtime_error("position " + to_string(pos) + " is in no sequence");
    return i;
}

}

GenomeConversion convertGenome(const Alignment *alignment, const string &name, const string &genomesDir,
                               const ConvertOptions &opts) {
    GenomeConversion out;
    const Genome *genome = alignment->openGenomeCheck(name);
    const Genome *parent = genome->getParent();
    vector<SeqInfo> seqs = sequencesOf(genome);
    vector<SeqInfo> pseqs;
    if (parent) pseqs = sequencesOf(parent);
    vector<api::SequenceDef> defs, pdefs;
    for (const SeqInfo &x : seqs) defs.push_back(api::SequenceDef{x.name, (int64_t)x.length});
    for (const SeqInfo &x : pseqs) pdefs.push_back(api::SequenceDef{x.name, (int64_t)x.length});
    api::GenomeWriterOptions wo;
    wo.chunkRuns = opts.chunkRuns;
    wo.chunkParentSpan = opts.chunkParentSpan;
    /* HAL's top segments come in genome order, one parent each at most: exactly the
     * writer's presorted promise, so runs stream to disk as they are merged */
    wo.presorted = true;
    wo.source = "halToFormat3";
    wo.checksums = false;   /* hal convert takes them once every genome is written */
    api::GenomeWriter w(genomesDir, name, parent ? parent->getName() : "", defs, pdefs, wo);
    for (const auto &kv : genome->getMetaData()->getMap()) w.setMeta(kv.first, kv.second);

    out.topSegments = genome->getNumTopSegments();
    if (parent && out.topSegments > 0) {
        TopSegmentIteratorPtr ti = genome->getTopSegmentIterator(0);
        BottomSegmentIteratorPtr bi = parent->getBottomSegmentIterator(0);
        size_t si = 0, pi = 0;
        for (hal_size_t k = 0; k < out.topSegments; k++, ti->toRight()) {
            const TopSegment *ts = ti->tseg();
            hal_index_t start = ts->getStartPosition();
            while (start >= seqs[si].start + (hal_index_t)seqs[si].length) si++;   /* segments never span sequences */
            if (!ts->hasParent()) continue;
            bi->toParent(ti);
            hal_index_t pstart = bi->bseg()->getStartPosition();
            hal_size_t len = ts->getLength();
            if (bi->bseg()->getLength() != len)
                throw runtime_error(name + ": top segment " + to_string(k) + " and its parent differ in length");
            if (pstart < pseqs[pi].start || pstart >= pseqs[pi].start + (hal_index_t)pseqs[pi].length)
                pi = sequenceAt(pseqs, pstart);
            w.addRun((uint32_t)si, start - seqs[si].start, (uint32_t)pi, pstart - pseqs[pi].start, (int64_t)len,
                     ts->getParentReversed());
            out.alignedSegments++;
            out.alignedBases += len;
        }
    }

    /* ---- the .2bit */
    if (opts.dna && genome->containsDNAArray()) {
        string buf;
        SequenceIteratorPtr sit = genome->getSequenceIterator();
        for (size_t i = 0; i < seqs.size(); i++, sit->toNext()) {
            const Sequence *s = sit->getSequence();
            w.beginDna();
            for (hal_size_t off = 0; off < seqs[i].length; off += 1 << 20) {
                hal_size_t n = min<hal_size_t>(1 << 20, seqs[i].length - off);
                s->getSubString(buf, off, n);
                w.appendDna(buf.data(), buf.size());
            }
            w.endDna();
        }
        out.wroteDna = true;
    }
    out.runs = w.finish().runs;
    alignment->closeGenome(genome);
    return out;
}

uint64_t checkGenome(const Alignment *alignment, const string &name, const Store &store, uint64_t nSamples,
                     uint64_t seed, ostream &log) {
    uint64_t bad = 0;
    auto complain = [&](const string &what) {
        if (bad++ < 10) log << "  " << name << ": " << what << "\n";
    };
    const Genome *genome = alignment->openGenomeCheck(name);
    const Genome *parent = genome->getParent();
    vector<SeqInfo> seqs = sequencesOf(genome), pseqs;
    if (parent) pseqs = sequencesOf(parent);
    string stem = "genomes/" + genomeFileStem(name);

    unique_ptr<lmap_file, void (*)(lmap_file *)> f(nullptr, lmap_close);
    if (store.has(stem + ".hgn")) f.reset(lmap_open_io(store.openMember(stem + ".hgn"), 1));
    if (!f) {
        complain("cannot open " + stem + ".hgn");
        return bad;
    }
    /* sequence tables, both axes */
    if (lmap_n_members(f.get(), 0) != seqs.size()) complain("axis a has the wrong number of sequences");
    else
        for (size_t i = 0; i < seqs.size(); i++) {
            const lmap_member *m = lmap_member_at(f.get(), 0, (uint32_t)i);
            if (seqs[i].name != m->name || seqs[i].length != m->length) complain("sequence " + seqs[i].name + " differs");
        }
    if (lmap_n_members(f.get(), 1) != pseqs.size()) complain("axis b has the wrong number of sequences");
    else
        for (size_t i = 0; i < pseqs.size(); i++) {
            const lmap_member *m = lmap_member_at(f.get(), 1, (uint32_t)i);
            if (pseqs[i].name != m->name || pseqs[i].length != m->length)
                complain("parent sequence " + pseqs[i].name + " differs");
        }
    unique_ptr<tb_file, void (*)(tb_file *)> dna(NULL, tb_close);
    if (genome->containsDNAArray()) {
        if (store.has(stem + ".2bit")) {
            lmap_io *io = store.openMember(stem + ".2bit");
            dna.reset(tb_open_io(io, 1));
            if (!dna) lmap_io_close(io);
        }
        if (!dna) complain("cannot open " + stem + ".2bit");
        else if (tb_n_seqs(dna.get()) != seqs.size()) complain(".2bit has the wrong number of sequences");
        else
            for (size_t i = 0; i < seqs.size(); i++)
                if (seqs[i].name != tb_name(dna.get(), (uint32_t)i) || tb_length(dna.get(), (uint32_t)i) != (int64_t)seqs[i].length)
                    complain(".2bit sequence " + seqs[i].name + " differs");
    }
    if (bad || genome->getSequenceLength() == 0) {
        alignment->closeGenome(genome);
        return bad;
    }

    /* random bases: parent position and strand, and DNA */
    uint64_t x = seed ? seed : 0x9E3779B97F4A7C15ull;
    auto rnd = [&x]() { x ^= x << 13; x ^= x >> 7; x ^= x << 17; return x; };
    hal_size_t glen = genome->getSequenceLength();
    TopSegmentIteratorPtr ti = parent && genome->getNumTopSegments() ? genome->getTopSegmentIterator(0) : TopSegmentIteratorPtr();
    BottomSegmentIteratorPtr bi = ti ? parent->getBottomSegmentIterator(0) : BottomSegmentIteratorPtr();
    string want;
    vector<char> got;
    for (uint64_t q = 0; q < nSamples; q++) {
        hal_index_t pos = (hal_index_t)(rnd() % glen);
        size_t si = sequenceAt(seqs, pos);
        int64_t local = pos - seqs[si].start;
        /* what HAL says */
        bool mapped = false, rev = false;
        size_t psi = 0;
        int64_t plocal = 0;
        if (ti) {
            ti->toSite(pos, false);
            const TopSegment *ts = ti->tseg();
            if (ts->hasParent()) {
                bi->toParent(ti);
                hal_index_t off = pos - ts->getStartPosition(), pstart = bi->bseg()->getStartPosition();
                rev = ts->getParentReversed();
                hal_index_t ppos = rev ? pstart + (hal_index_t)ts->getLength() - 1 - off : pstart + off;
                psi = sequenceAt(pseqs, ppos);
                plocal = ppos - pseqs[psi].start;
                mapped = true;
            }
        }
        /* what the .hgn says */
        lmap_hit *hits = NULL;
        size_t n = 0;
        if (lmap_query_a(f.get(), (uint32_t)si, local, local + 1, &hits, &n, NULL) != 0) {
            complain("query failed at " + seqs[si].name + ":" + to_string(local));
            continue;
        }
        string where = seqs[si].name + ":" + to_string(local);
        if (!mapped && n != 0) complain(where + " is unaligned in HAL but mapped in the .hgn");
        else if (mapped && n != 1) complain(where + " has " + to_string(n) + " parent positions in the .hgn");
        else if (mapped) {
            const lmap_hit &h = hits[0];
            int64_t bpos = h.strand ? h.b + h.len - 1 - (local - h.a) : h.b + (local - h.a);
            if (h.b_member != psi || bpos != plocal || (h.strand != 0) != rev)
                complain(where + ": HAL gives " + pseqs[psi].name + ":" + to_string(plocal) + (rev ? " -" : " +") +
                         ", the .hgn " + (h.b_member < pseqs.size() ? pseqs[h.b_member].name : "#" + to_string(h.b_member)) + ":" + to_string(bpos) + (h.strand ? " -" : " +"));
        }
        free(hits);
        /* DNA, a short window */
        if (dna) {
            int64_t w = min<int64_t>(100, (int64_t)seqs[si].length - local);
            genome->getSequenceBySite(pos)->getSubString(want, (hal_size_t)local, (hal_size_t)w);
            got.resize((size_t)w);
            if (tb_fetch(dna.get(), (uint32_t)si, local, local + w, got.data(), 0) != 0 ||
                string(got.begin(), got.end()) != want)
                complain("DNA differs at " + where);
        }
    }
    alignment->closeGenome(genome);
    return bad;
}
}
