/*
 * HAL format 3: the C interface.  See hal3/hal3_c.h.
 */
#include "hal3/hal3_c.h"
#include "hal3/hal3Cactus.h"

#include "hal3Columns.h"
#include "hal3Store.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct hal3_c2h {
    unique_ptr<hal3::api::C2hImporter> importer;
};

struct hal3_columns {
    unique_ptr<hal3::Alignment> aln;
    unique_ptr<hal3::UniversalColumns> cols;   /* (declared after aln: it goes first) */
    vector<hal3::ColumnChunk> chunks;          /* per worker */
};

namespace {

void setError(char **err, const char *what) {
    if (err) *err = strdup(what ? what : "unknown error");
}

/* run f, turning an exception into -1 and *err */
template <class F>
int guarded(char **err, F f) {
    try {
        f();
        return 0;
    } catch (exception &e) {
        setError(err, e.what());
    } catch (...) {
        setError(err, "unknown error");
    }
    return -1;
}
}

extern "C" {

hal3_c2h *hal3_c2h_open(const char *fragmentDir, const char *newick, int presorted, char **err) {
    hal3_c2h *h = nullptr;
    int rc = guarded(err, [&]() {
        hal3::api::GenomeWriterOptions o;
        o.presorted = presorted != 0;
        o.source = "cactus";
        unique_ptr<hal3_c2h> p(new hal3_c2h());
        p->importer.reset(new hal3::api::C2hImporter(fragmentDir, newick, o));
        h = p.release();
    });
    return rc == 0 ? h : nullptr;
}

int hal3_c2h_declare(hal3_c2h *h, const char *genome, const char *name, int ancestor, int64_t length, char **err) {
    return guarded(err, [&]() { h->importer->declare(genome, name, ancestor != 0, length); });
}

int hal3_c2h_sequence(hal3_c2h *h, const char *genome, const char *name, const char *segments, size_t segmentsLength,
                      const char *dna, int64_t dnaLength, char **err) {
    return guarded(err, [&]() { h->importer->sequence(genome, name, segments, segmentsLength, dna, dnaLength); });
}

int hal3_c2h_finish(hal3_c2h *h, char **err) {
    int rc = guarded(err, [&]() { h->importer->finish(); });
    delete h;
    return rc;
}

void hal3_c2h_free(hal3_c2h *h) { delete h; }

int hal3_pack(const char *in, const char *out, char **err) {
    return guarded(err, [&]() { hal3::packAlignment(in, out); });
}

void hal3_columns_defaults(hal3_columns_options *opts) {
    memset(opts, 0, sizeof *opts);
    opts->length = -1;
    opts->ancestorBases = 1;
    opts->chunk = hal3::ColumnOptions().chunk;
    opts->workers = 1;
}

hal3_columns *hal3_columns_open(const char *alignment, const hal3_columns_options *o, char **err) {
    hal3_columns *c = nullptr;
    int rc = guarded(err, [&]() {
        hal3_columns_options d;
        hal3_columns_defaults(&d);
        if (!o) o = &d;
        if (o->workers < 1 || o->chunk < 1) throw runtime_error("workers and chunk must be positive");
        unique_ptr<hal3_columns> p(new hal3_columns());
        p->aln.reset(new hal3::Alignment(alignment));
        hal3::ColumnOptions co;
        if (o->rootGenome) co.rootGenome = o->rootGenome;
        if (o->excludeGenomes) {
            string s = o->excludeGenomes;
            for (size_t a = 0; a <= s.size();) {
                size_t b = s.find(',', a);
                if (b == string::npos) b = s.size();
                if (b > a) co.exclude.push_back(s.substr(a, b - a));
                a = b + 1;
            }
        }
        co.leafInsertions = o->leafInsertions != 0;
        co.ancestorBases = o->ancestorBases != 0;
        co.chunk = o->chunk;
        vector<hal3::RefInterval> ranges;
        if (o->refSequence) {
            const hal3::Genome *b = co.rootGenome.empty() ? p->aln->root() : p->aln->genome(co.rootGenome);
            if (!b) throw runtime_error("no genome " + co.rootGenome);
            int64_t i = b->findSequence(o->refSequence);
            if (i < 0) throw runtime_error(string("no sequence ") + o->refSequence + " in " + b->name());
            const int64_t len = b->sequenceLength((uint32_t)i);
            if (o->start < 0 || o->start >= len) throw runtime_error("start is outside " + string(o->refSequence));
            const int64_t e = o->length < 0 ? len : min(len, o->start + o->length);
            if (o->start >= e) throw runtime_error("the range is empty");
            ranges.push_back(hal3::RefInterval{(uint32_t)i, o->start, e});
        } else if (o->start != 0 || o->length >= 0) {
            throw runtime_error("start and length need refSequence");
        }
        p->cols.reset(new hal3::UniversalColumns(*p->aln, ranges, co, o->workers));
        p->chunks.resize((size_t)o->workers);
        c = p.release();
    });
    return rc == 0 ? c : nullptr;
}

size_t hal3_columns_chunks(const hal3_columns *c) { return c->cols->chunks(); }
int64_t hal3_columns_columns(const hal3_columns *c) { return c->cols->columns(); }
uint32_t hal3_columns_genomes(const hal3_columns *c) { return c->cols->genomes(); }
const char *hal3_columns_genome_name(const hal3_columns *c, uint32_t i) { return c->cols->genome(i)->name().c_str(); }
int32_t hal3_columns_genome_parent(const hal3_columns *c, uint32_t i) { return c->cols->parentOf(i); }
int hal3_columns_genome_is_leaf(const hal3_columns *c, uint32_t i) { return c->cols->genome(i)->children().empty() ? 1 : 0; }

double hal3_columns_branch_length(const hal3_columns *c, uint32_t i) {
    const hal3::GenomeEntry *e = c->aln->manifest().find(c->cols->genome(i)->name());
    return e ? e->branchLength : -1;
}

int hal3_columns_chunk(hal3_columns *c, size_t k, int w, hal3_column_chunk *out, char **err) {
    return guarded(err, [&]() {
        if (w < 0 || w >= (int)c->chunks.size()) throw runtime_error("no worker " + to_string(w));
        hal3::ColumnChunk &ch = c->chunks[(size_t)w];
        c->cols->chunk(k, w, ch);
        out->axisName = ch.axisName->c_str();
        out->axisStart = ch.axisStart;
        out->cols = ch.cols;
        out->colStart = ch.colStart.data();
        out->genome = ch.genome.data();
        out->parent = ch.parent.data();
        out->piece = ch.piece.data();
        out->base = ch.base.data();
        out->entries = (int64_t)ch.base.size();
    });
}

void hal3_columns_close(hal3_columns *c) { delete c; }
}
