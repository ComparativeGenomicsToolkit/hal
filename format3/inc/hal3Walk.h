/*
 * HAL format 3: the tree walk behind reference-anchored output (hal maf, hal coverage).
 *
 * Lifts intervals of a reference genome to every genome in scope at once.
 *
 * Default: up the reference's lineage to the scope top, collecting at each genome the
 * bases that go no higher ("tops"), then one pass down from the top that descends from
 * every top into its whole subtree.  That is hal2maf's column: every descendant of the
 * highest ancestral base a reference base reaches, paralogs included.
 *
 * onlyOrthologs: at each ancestor on the reference's lineage, descend only into the
 * children off the lineage (and into the reference's own subtree): liftover to each
 * genome's common ancestor with the reference, as halLiftover and halCoverage.
 *
 * dupes off (--noDupes): up only from canonical copies, down only to them.
 *
 * Scope: the spanning tree of the printed genomes -- the targets, or the subtree of
 * rootGenome, or every genome -- and the reference.  Only printed genomes get pieces;
 * ancestors are left out with ancestors = false, excluded genomes always, and the
 * reference's own other copies with refSingleCopy.  Leaving a genome out of the printed
 * set changes no other genome's pieces: the walk still goes through it.
 *
 * Also here, for the commands built on the walk: cutting reference intervals into
 * chunks, and running chunks on several threads with results taken in order.
 */
#ifndef _HAL3WALK_H
#define _HAL3WALK_H

#include "hal3Alignment.h"
#include "hal3Liftover.h"

#include <algorithm>
#include <condition_variable>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace hal3 {

struct WalkOptions {
    bool ancestors = true;
    bool dupes = true;
    bool onlyOrthologs = false;
    std::vector<std::string> targets;   /* empty: every genome */
    std::string rootGenome;
    std::vector<std::string> exclude;   /* never printed */
    bool refSingleCopy = false;         /* the reference's other copies not printed */
};

class TreeWalk {
  public:
    /* Throws std::runtime_error on unknown genomes or bad combinations. */
    TreeWalk(const Alignment &alignment, const Genome *ref, const WalkOptions &opts);

    /* Lift intervals of the reference sequence refSeq (sorted, disjoint) to every
     * printed genome.  pieces[g]: source = the reference, target = g.  The reference's
     * own identity pieces are left out; its paralogous copies are not. */
    void lift(uint32_t refSeq, const std::vector<std::pair<int64_t, int64_t>> &intervals,
              std::map<const Genome *, std::vector<LiftPiece>> &pieces);

    const Genome *ref() const { return _ref; }
    const Genome *top() const { return _top; }
    /* the reference, its parent, ... up to top() */
    const std::vector<const Genome *> &refPath() const { return _refPath; }
    const std::set<const Genome *> &printed() const { return _printed; }
    const WalkOptions &options() const { return _opts; }
    CursorCache &cursors() { return _cursors; }

  private:
    const Alignment &_aln;
    const Genome *_ref;
    WalkOptions _opts;
    std::set<const Genome *> _scope, _printed;
    std::vector<const Genome *> _refPath;
    const Genome *_top = nullptr;
    CursorCache _cursors;

    void descend(const Genome *v, std::vector<LiftPiece> &in, std::map<const Genome *, std::vector<LiftPiece>> &tops,
                 std::map<const Genome *, std::vector<LiftPiece>> &pieces);
};

/* Merge pieces contiguous in both source and target. */
void coalescePieces(std::vector<LiftPiece> &v);

/* [start, end) of reference sequence seq. */
struct RefInterval {
    uint32_t seq;
    int64_t start, end;
};

/* Each interval cut every `size` bases from its start, in order. */
std::vector<RefInterval> chunkIntervals(const std::vector<RefInterval> &intervals, int64_t size);

/* produce(k) for k = 0 .. n-1 on `threads` threads, and consume(k, result) on the
 * calling thread in order of k.  A thread takes the next k only while fewer than
 * 2 x threads results are taken and not yet consumed, which bounds what waits.  The
 * first exception, from either, stops the rest and is rethrown once every thread has
 * finished.  With one thread, each k is produced and consumed in turn on the calling
 * thread. */
template <class R>
void orderedParallel(size_t n, int threads, const std::function<R(size_t)> &produce,
                     const std::function<void(size_t, R &)> &consume) {
    if (threads <= 1 || n <= 1) {
        for (size_t k = 0; k < n; k++) {
            R r = produce(k);
            consume(k, r);
        }
        return;
    }
    std::mutex m;
    std::condition_variable cv;
    std::vector<std::unique_ptr<R>> done(n);
    size_t taken = 0, consumed = 0, ahead = 2 * (size_t)threads;
    bool stop = false;
    std::exception_ptr error;
    auto fail = [&]() {   /* with m held */
        if (!error) error = std::current_exception();
        stop = true;
        cv.notify_all();
    };
    auto work = [&]() {
        for (;;) {
            size_t k;
            {
                std::unique_lock<std::mutex> lock(m);
                cv.wait(lock, [&]() { return stop || taken >= n || taken < consumed + ahead; });
                if (stop || taken >= n) return;
                k = taken++;
            }
            try {
                std::unique_ptr<R> r(new R(produce(k)));
                std::lock_guard<std::mutex> lock(m);
                done[k] = std::move(r);
                cv.notify_all();
            } catch (...) {
                std::lock_guard<std::mutex> lock(m);
                fail();
                return;
            }
        }
    };
    std::vector<std::thread> pool;
    try {
        for (size_t t = 0; t < std::min<size_t>((size_t)threads, n); t++) pool.emplace_back(work);
    } catch (...) {   /* a thread could not be started: stop the others, then report it */
        {
            std::lock_guard<std::mutex> lock(m);
            stop = true;
            cv.notify_all();
        }
        for (std::thread &t : pool) t.join();
        throw;
    }
    for (size_t k = 0; k < n; k++) {
        std::unique_ptr<R> r;
        {
            std::unique_lock<std::mutex> lock(m);
            cv.wait(lock, [&]() { return stop || done[k]; });
            if (stop) break;
            r = std::move(done[k]);
        }
        try {
            consume(k, *r);
        } catch (...) {
            std::lock_guard<std::mutex> lock(m);
            fail();
            break;
        }
        std::lock_guard<std::mutex> lock(m);
        consumed = k + 1;
        cv.notify_all();
    }
    for (std::thread &t : pool) t.join();
    if (error) std::rethrow_exception(error);
}
}

#endif
