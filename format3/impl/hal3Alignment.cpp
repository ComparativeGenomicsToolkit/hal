/*
 * HAL format 3: reading an alignment.  See hal3Alignment.h.
 */
#include "hal3Alignment.h"

#include <stdexcept>

using namespace std;

namespace hal3 {

/* Nothing is thrown out of call_once: an exception leaving it relies on pthread_once
 * unwinding, which some runtimes (musl, ThreadSanitizer) do not support -- the other
 * threads would wait forever.  A failure is recorded and thrown after. */
void Genome::open() const {
    call_once(_opened, [this]() {
        string path = _alignment->path() + ":" + _entry->hgnPath;
        try {
            lmap_io *io = _alignment->store().openMember(_entry->hgnPath);
            lmap_file *f = lmap_open_io(io, 1);
            if (!f) throw runtime_error("cannot open genome file " + path);
            const char *g = lmap_meta_get(f, "hal.genome");
            if (!g || _entry->name != g) {
                lmap_close(f);
                throw runtime_error(path + " is not the genome file of " + _entry->name);
            }
            uint32_t n = lmap_n_members(f, 0);
            _starts.resize(n + 1);
            _starts[0] = 0;
            for (uint32_t i = 0; i < n; i++) _starts[i + 1] = _starts[i] + (int64_t)lmap_member_at(f, 0, i)->length;
            _hgn = f;
        } catch (exception &e) {
            _openError = e.what();
        }
    });
    if (!_hgn) throw runtime_error(_openError.empty() ? "genome file of " + _entry->name + " failed to open" : _openError);
}

lmap_file *Genome::hgn() const {
    open();
    return _hgn;
}

uint32_t Genome::numSequences() const {
    open();
    return lmap_n_members(_hgn, 0);
}

const char *Genome::sequenceName(uint32_t i) const {
    open();
    return lmap_member_at(_hgn, 0, i)->name;
}

int64_t Genome::sequenceLength(uint32_t i) const {
    open();
    return _starts[i + 1] - _starts[i];
}

int64_t Genome::sequenceStart(uint32_t i) const {
    open();
    return _starts[i];
}

int64_t Genome::findSequence(const string &name) const {
    open();
    return lmap_member_by_name(_hgn, 0, name.c_str());
}

Alignment::Alignment(const string &path) : _path(path) {
    _store = Store::open(path);
    _manifest = Manifest::parse(_store->read(MANIFEST_NAME), path);
    for (const GenomeEntry &e : _manifest.genomes) {
        unique_ptr<Genome> g(new Genome());
        g->_alignment = this;
        g->_entry = &e;
        _byName[e.name] = g.get();
        _genomes.push_back(move(g));
    }
    /* the roster lists parents before children (hal3Validate checks it) */
    for (auto &g : _genomes) {
        if (g->_entry->parent.empty()) continue;
        auto p = _byName.find(g->_entry->parent);
        if (p == _byName.end()) throw runtime_error("genome " + g->name() + " names unknown parent " + g->_entry->parent);
        g->_parent = p->second;
        g->_depth = p->second->_depth + 1;
        const_cast<Genome *>(p->second)->_children.push_back(g.get());
    }
    if (_genomes.empty() || _genomes[0]->_parent) throw runtime_error(_path + ": the root must be listed first");
}

Alignment::~Alignment() {
    for (auto &g : _genomes) {
        if (g->_hgn) lmap_close(g->_hgn);
        if (g->_dna) lmap_io_close(g->_dna);
    }
}

tb_file *Alignment::openDna(const Genome *g, bool shared) const {
    const GenomeEntry *e = _manifest.find(g->name());
    if (!e || e->dnaPath.empty()) throw runtime_error("genome " + g->name() + " has no DNA in this alignment");
    lmap_io *io = nullptr;
    if (shared) {
        call_once(g->_dnaOpened, [&]() {   /* nothing thrown out of call_once: see Genome::open */
            try {
                g->_dna = _store->openMember(e->dnaPath);
            } catch (exception &x) {
                g->_dnaError = x.what();
            }
        });
        if (!g->_dna) throw runtime_error(g->_dnaError.empty() ? "cannot read " + e->dnaPath + " of " + _path : g->_dnaError);
        io = lmap_io_slice(g->_dna, 0, lmap_io_size(g->_dna));
    } else {
        io = _store->openMember(e->dnaPath);
    }
    tb_file *f = io ? tb_open_io(io, 1) : nullptr;
    if (!f) {
        if (io) lmap_io_close(io);
        throw runtime_error("cannot read " + e->dnaPath + " of " + _path);
    }
    return f;
}

const Genome *Alignment::genome(const string &name) const {
    auto i = _byName.find(name);
    return i == _byName.end() ? nullptr : i->second;
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
}
