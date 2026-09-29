/*
 * HAL format 3: the C interface.  See hal3/hal3_c.h.
 */
#include "hal3/hal3_c.h"
#include "hal3/hal3Cactus.h"

#include "hal3Store.h"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <string>

using namespace std;

struct hal3_c2h {
    unique_ptr<hal3::api::C2hImporter> importer;
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
}
