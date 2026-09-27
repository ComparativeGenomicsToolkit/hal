/*
 * hal likelihoods, hal posteriors: per-base likelihoods of the ancestral genomes (their
 * "lik" members, hal3Likelihood.h), a reference implementation.
 *
 * Model: Jukes-Cantor (JC69), the manifest's branch lengths (below 1e-6 counted as 1e-6,
 * so that conflicting children never make a vector all zero), a uniform prior.  A leaf's
 * base is certain (N: no information); every aligned copy of a base in a child counts as
 * independent evidence (paralogy); a base with no copy in a child learns nothing from it.
 *
 * hal likelihoods -- bottom up, children before parents (Felsenstein's pruning): for
 * each base a of an ancestor, over every child c and every copy of a in c,
 *   L(a) = prod  sum_b P_c(b | a) L_c(b),   normalized,
 * L_c being a leaf child's base or an ancestral child's own L as just written (so each
 * level uses its children's stored values, as building an alignment bottom-up would).
 *
 * hal posteriors -- top down, parents before children: a child base b whose parent base
 * is p gets, from the parent's posterior P_p with this copy's own message m_c divided
 * out,
 *   U(b) = sum_a P_c(b | a) P_p(a) / m_c(a),    posterior(b) ~ L(b) U(b)
 * (U uniform for a base with no parent), and the root's posterior is its L.  Each
 * ancestral base is then set to the most probable (unless --keepBases), keeping its
 * case, and the posteriors replace the likelihoods.
 *
 * Genomes are done one after another, each in windows on several threads (the chunks
 * encoded there too, and written in order, so the files are the same for any number of
 * threads), into scratch files beside the output; then one update (hal3Update.h)
 * attaches them: a new alignment, or --inPlace one crash-safe change.  Everything that
 * could stop the update is checked before the work starts, and the input must not have
 * changed by the end.
 */
#include "hal3Alignment.h"
#include "hal3Cli.h"
#include "hal3Likelihood.h"
#include "hal3Liftover.h"
#include "hal3Update.h"
#include "hal3Walk.h"
#include "twobit64.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <map>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

namespace hal3 {

namespace {

const uint32_t CHUNK = 65536;                  /* bases per .hlk chunk */
const int64_t WINDOW = (int64_t)16 * CHUNK;    /* a thread's unit: whole chunks */
const uint64_t SAMPLE = (uint64_t)1 << 22;     /* bases whose codes choose the dictionary */
const double MIN_BRANCH = 1e-6;
const int64_t SPAN_GAP = 4096;                 /* pieces this close are read as one span */
const int64_t SPAN_MAX = (int64_t)1 << 18;

/* JC69 over a branch of length t: P(b | a) is `same` for b = a, `other` otherwise. */
struct Jc {
    double same, other;
    explicit Jc(double t) {
        double e = exp(-4.0 * max(t, MIN_BRANCH) / 3.0);
        same = 0.25 + 0.75 * e;
        other = 0.25 - 0.25 * e;
    }
};

struct Options {
    int threads = 1, q = 10;
    bool inPlace = false, directory = false, keepBases = false;
    string in, out;
};

/* Scratch files beside the output (or the input, in place), removed at the end. */
struct Scratch {
    string dir;
    vector<string> files;
    explicit Scratch(const string &near) : dir(near + ".lik-tmp") {
        if (mkdir(dir.c_str(), 0777) != 0)
            throw runtime_error("cannot create " + dir + (errno == EEXIST ? " (left by an interrupted run? remove it)" : ""));
    }
    ~Scratch() {
        for (const string &f : files) {
            unlink(f.c_str());
            unlink((f + ".partial").c_str());
        }
        rmdir(dir.c_str());
    }
    string file(const string &name) {
        files.push_back(dir + "/" + name);
        return files.back();
    }
};

/* One genome's DNA and, if it has them, likelihoods: what a worker reads of it. */
struct GenomeData {
    string dnaFile;                 /* a local .2bit replacing the alignment's, or "" */
    unique_ptr<LikFile> lik;        /* its (new or stored) likelihoods, or none */
};

/* Per window: cursors, DNA handles, a likelihood chunk cache per genome. */
struct Worker {
    const Alignment &aln;
    const map<const Genome *, GenomeData> &data;
    CursorCache cursors;
    map<const Genome *, tb_file *> dna;
    map<const Genome *, LikCache> caches;
    string bases;
    Worker(const Alignment &a, const map<const Genome *, GenomeData> &d) : aln(a), data(d), cursors((size_t)1 << 20, 256) {}
    ~Worker() {
        for (auto &kv : dna) tb_close(kv.second);
    }
    tb_file *dnaOf(const Genome *g) {
        auto it = dna.find(g);
        if (it != dna.end()) return it->second;
        tb_file *f;
        auto d = data.find(g);
        if (d != data.end() && !d->second.dnaFile.empty()) {
            f = tb_open(d->second.dnaFile.c_str());
            if (!f) throw runtime_error("cannot read " + d->second.dnaFile);
        } else {
            f = aln.openDna(g, true);
        }
        return dna[g] = f;
    }
    /* bases [s, e) of g's sequence seq, forward (a format 3 .2bit holds the sequence
     * table's sequences in its order) */
    const string &fetch(const Genome *g, uint32_t seq, int64_t s, int64_t e) {
        bases.resize((size_t)(e - s));
        if (e > s && tb_fetch(dnaOf(g), seq, s, e, &bases[0], 0) != 0)
            throw runtime_error("cannot read DNA of " + g->name() + "." + g->sequenceName(seq));
        return bases;
    }
    /* g's vectors for [s, e) of seq, forward, 4 per base: its likelihoods, or its bases
     * as certain (N: 1/4 each) */
    void vectors(const Genome *g, uint32_t seq, int64_t s, int64_t e, vector<double> &out) {
        out.resize((size_t)(4 * (e - s)));
        const string &b = fetch(g, seq, s, e);
        auto d = data.find(g);
        if (d != data.end() && d->second.lik) {
            d->second.lik->probabilities(seq, s, e, b.data(), false, out.data(), &caches[g]);
            return;
        }
        for (size_t i = 0; i < b.size(); i++) {
            int x = likBase(b[i]);
            double *v = &out[4 * i];
            for (int k = 0; k < 4; k++) v[k] = x == 4 ? 0.25 : (k == x ? 1.0 : 0.0);
        }
    }
    /* each(piece, v) for every piece, v being g's vectors for [tgtStart, tgtStart +
     * length) of its tgtSeq: the pieces are taken in order along g and read a span at a
     * time where they lie close together, so a chunk or a stretch of DNA is read once
     * however many pieces fall in it (many, where the bases have many copies) */
    template <class F>
    void forPieces(const Genome *g, const vector<LiftPiece> &pieces, F each) {
        vector<size_t> order(pieces.size());
        for (size_t i = 0; i < order.size(); i++) order[i] = i;
        sort(order.begin(), order.end(), [&](size_t x, size_t y) {
            return pieces[x].tgtSeq != pieces[y].tgtSeq ? pieces[x].tgtSeq < pieces[y].tgtSeq : pieces[x].tgtStart < pieces[y].tgtStart;
        });
        vector<double> v;
        for (size_t i = 0; i < order.size();) {
            const LiftPiece &a = pieces[order[i]];
            int64_t s = a.tgtStart, e = a.tgtStart + a.length;
            size_t j = i + 1;
            for (; j < order.size(); j++) {
                const LiftPiece &b = pieces[order[j]];
                if (b.tgtSeq != a.tgtSeq || b.tgtStart > e + SPAN_GAP || max(e, b.tgtStart + b.length) - s > SPAN_MAX) break;
                e = max(e, b.tgtStart + b.length);
            }
            vectors(g, a.tgtSeq, s, e, v);
            for (size_t k = i; k < j; k++) {
                const LiftPiece &p = pieces[order[k]];
                each(p, &v[(size_t)(4 * (p.tgtStart - s))]);
            }
            i = j;
        }
    }
};

vector<pair<string, int64_t>> sequences(const Genome *g) {
    vector<pair<string, int64_t>> v;
    for (uint32_t i = 0; i < g->numSequences(); i++) v.push_back({g->sequenceName(i), g->sequenceLength(i)});
    return v;
}

bool sameSequences(const LikFile &f, const Genome *g) {
    if (f.numSequences() != g->numSequences()) return false;
    for (uint32_t i = 0; i < f.numSequences(); i++)
        if (f.sequenceName(i) != g->sequenceName(i) || f.sequenceLength(i) != g->sequenceLength(i)) return false;
    return true;
}

void normalize(double *v) {
    double s = v[0] + v[1] + v[2] + v[3];
    if (!(s > 0) || !isfinite(s)) {
        v[0] = v[1] = v[2] = v[3] = 0.25;
        return;
    }
    for (int k = 0; k < 4; k++) v[k] /= s;
}

/* A window's result: its codes, and (top down) its new bases. */
struct WindowOut {
    vector<LikCode> codes;
    string bases;
};

/* One genome's pass: compute(k) gives window k's codes (and bases); the dictionary comes
 * from the first SAMPLE bases' windows; every window's chunks are encoded on the
 * threads and written in order, and its bases, if any, to tw. */
void runGenome(const Genome *g, const vector<RefInterval> &wins, const string &file, const map<string, string> &meta,
               const Options &o, const function<WindowOut(size_t)> &compute, tb_writer *tw, const string &dnaFile) {
    size_t ns = 0;
    for (uint64_t nb = 0; ns < wins.size() && nb < SAMPLE; ns++) nb += (uint64_t)(wins[ns].end - wins[ns].start);
    vector<WindowOut> sample(ns);
    orderedParallel<WindowOut>(ns, o.threads, compute, [&](size_t k, WindowOut &r) { sample[k] = move(r); });
    vector<LikCode> all;
    for (const WindowOut &w : sample) all.insert(all.end(), w.codes.begin(), w.codes.end());
    LikDictionary dict = LikDictionary::choose(all.data(), all.size());
    vector<LikCode>().swap(all);
    LikWriter writer(file, sequences(g), meta, dict, o.q, CHUNK);
    struct Encoded {
        vector<LikChunk> chunks;
        string bases;
    };
    auto check = [&](int rc) {
        if (rc != 0) throw runtime_error(dnaFile + ": " + (tb_writer_error(tw) ? tb_writer_error(tw) : "write failed"));
    };
    int64_t nextSeq = 0;   /* the .2bit: the sequence to begin next */
    auto emptyUpTo = [&](int64_t s) {   /* sequences with no windows (no bases) up to s */
        for (; nextSeq < s; nextSeq++) {
            check(tb_writer_begin(tw));
            check(tb_writer_end(tw));
        }
    };
    int64_t openSeq = -1;
    orderedParallel<Encoded>(
        wins.size(), o.threads,
        [&](size_t k) {
            WindowOut r = k < ns ? move(sample[k]) : compute(k);
            Encoded e;
            for (size_t at = 0; at < r.codes.size(); at += CHUNK)
                e.chunks.push_back(encodeLikChunk(r.codes.data() + at, min<size_t>(CHUNK, r.codes.size() - at), dict));
            e.bases = move(r.bases);
            return e;
        },
        [&](size_t k, Encoded &e) {
            for (const LikChunk &c : e.chunks) writer.addChunk(c);
            if (!tw) return;
            int64_t seq = wins[k].seq;
            if (seq != openSeq) {
                if (openSeq >= 0) check(tb_writer_end(tw));
                emptyUpTo(seq);
                check(tb_writer_begin(tw));
                openSeq = seq;
                nextSeq = seq + 1;
            }
            check(tb_writer_append(tw, e.bases.data(), e.bases.size()));
        });
    writer.close();
    if (tw) {
        if (openSeq >= 0) check(tb_writer_end(tw));
        emptyUpTo(g->numSequences());
    }
}

vector<RefInterval> windows(const Genome *g) {
    vector<RefInterval> all;
    for (uint32_t i = 0; i < g->numSequences(); i++) all.push_back(RefInterval{i, 0, g->sequenceLength(i)});
    return chunkIntervals(all, WINDOW);
}

vector<const Genome *> preorderAll(const Alignment &aln) {
    vector<const Genome *> out, stack{aln.root()};
    while (!stack.empty()) {
        const Genome *g = stack.back();
        stack.pop_back();
        out.push_back(g);
        for (auto c = g->children().rbegin(); c != g->children().rend(); ++c) stack.push_back(*c);
    }
    return out;
}

void parse(int argc, char **argv, const char *usage, Options &o, bool posteriors) {
    enum { THREADS = 1000, STEPS, IN_PLACE, DIRECTORY, KEEP_BASES };
    static const option longOpts[] = {{"threads", required_argument, nullptr, THREADS},
                                      {"stepsPerDecade", required_argument, nullptr, STEPS},
                                      {"inPlace", no_argument, nullptr, IN_PLACE},
                                      {"directory", no_argument, nullptr, DIRECTORY},
                                      {"keepBases", no_argument, nullptr, KEEP_BASES},
                                      {"help", no_argument, nullptr, 'h'},
                                      {nullptr, 0, nullptr, 0}};
    int c;
    int64_t threads = 1, q = 10;
    while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
        switch (c) {
        case THREADS: threads = parseInt("threads", optarg); break;
        case STEPS: q = parseInt("stepsPerDecade", optarg); break;
        case IN_PLACE: o.inPlace = true; break;
        case DIRECTORY: o.directory = true; break;
        case KEEP_BASES:
            if (!posteriors) throw CliError("bad option");
            o.keepBases = true;
            break;
        case 'h': cout << usage; exit(0);
        default: throw CliError("bad option");
        }
    }
    int want = o.inPlace ? 1 : 2;
    if (argc - optind != want) throw CliError("expected " + to_string(want) + " argument" + (want > 1 ? "s" : "") + ", got " + to_string(argc - optind));
    if (threads < 1 || threads > 4096) throw CliError("--threads must be 1 to 4096");
    if (q < 1 || q > LIK_MAX_STEPS) throw CliError("--stepsPerDecade must be 1 to " + to_string(LIK_MAX_STEPS));
    if (o.inPlace && o.directory) throw CliError("--inPlace changes the input in its own form; --directory is for new output");
    o.threads = (int)threads;
    o.q = (int)q;
    o.in = argv[optind];
    if (!o.inPlace) o.out = argv[optind + 1];
}

/* What would stop the final update, or the model, checked before any work. */
void preflight(const Options &o, const Alignment &aln) {
    struct stat st;
    string target = o.inPlace ? o.in : o.out;
    if (!o.inPlace && stat(o.out.c_str(), &st) == 0) throw runtime_error(o.out + " already exists");
    if (stat((target + ".update-tmp").c_str(), &st) == 0)
        throw runtime_error(target + ".update-tmp exists (left by an interrupted update? remove it)");
    if (o.inPlace) {
        if (access(o.in.c_str(), W_OK) != 0) throw runtime_error("cannot change " + o.in + ": " + strerror(errno));
    } else {
        size_t k = o.out.rfind('/');
        string dir = k == string::npos ? "." : k == 0 ? "/" : o.out.substr(0, k);
        if (access(dir.c_str(), W_OK) != 0) throw runtime_error("cannot write in " + dir + ": " + strerror(errno));
    }
    for (const GenomeEntry &e : aln.manifest().genomes) {
        if (!e.parent.empty() && !(e.branchLength >= 0))
            throw runtime_error("genome " + e.name + " has no branch length; the model needs one");
        if (e.dnaPath.empty()) throw runtime_error("genome " + e.name + " has no DNA; likelihoods need every genome's");
    }
}

/* Attach the files -- every ancestor's likelihoods (and new DNA), leaves' old ones
 * dropped -- in one update, if the input is as it was. */
void attach(const Options &o, const string &manifestBefore, const map<const Genome *, GenomeData> &data,
            const map<const Genome *, string> &likFiles, const char *name) {
    Update u(o.in);
    if (Store::open(o.in)->read(MANIFEST_NAME) != manifestBefore)
        throw runtime_error(o.in + " changed while " + name + " ran; nothing written");
    for (auto &kv : likFiles) {
        u.setMember(kv.first->name(), "lik", kv.second);
        auto d = data.find(kv.first);
        if (d != data.end() && !d->second.dnaFile.empty()) u.setDna(kv.first->name(), d->second.dnaFile);
    }
    Manifest m = Manifest::parse(manifestBefore, o.in);
    for (const GenomeEntry &e : m.genomes) {
        bool leaf = true;
        for (const GenomeEntry &c : m.genomes)
            if (c.parent == e.name) leaf = false;
        if (leaf && e.member("lik")) u.dropMember(e.name, "lik");
    }
    if (o.inPlace) u.writeInPlace();
    else u.write(o.out, o.directory);
    for (const string &l : u.log())
        if (l.compare(0, 4, "set ") != 0) cerr << "hal " << name << ": " << l << endl;
    cerr << "hal " << name << ": " << likFiles.size() << " ancestral genomes" << (o.inPlace ? ", in place" : "") << endl;
}
}

static const char *LIKELIHOODS_USAGE =
    "usage: hal likelihoods [options] <in.hal> <out.hal>\n"
    "       hal likelihoods --inPlace [options] <in.hal>\n"
    "For every base of every ancestral genome, the likelihood of the data below it given each\n"
    "base (A, C, G, T), normalized: Felsenstein's pruning under JC69 with the tree's branch\n"
    "lengths, every aligned copy in a child counting.  Stored as each ancestor's \"subtree\"\n"
    "likelihoods (its lik member), replacing any it has; leaves' likelihoods are dropped.  See\n"
    "hal posteriors for the top-down pass.\n"
    "options:\n"
    "  --threads N          windows of a genome to do at once (default 1)\n"
    "  --stepsPerDecade Q   precision of the stored values: log10 ratios to 1/Q, 1 to 25 (default\n"
    "                       10); a file spans 254/Q decades, smaller probabilities stored as that\n"
    "  --inPlace            change in.hal itself (one crash-safe transaction for a single file)\n"
    "  --directory          write out.hal as a directory\n"
    "  -h, --help\n";

int main_likelihoods(int argc, char **argv) {
    return runCommand("likelihoods", [&]() {
        Options o;
        parse(argc, argv, LIKELIHOODS_USAGE, o, false);
        Alignment aln(o.in);
        string manifestBefore = aln.store().read(MANIFEST_NAME);
        preflight(o, aln);
        Scratch scratch(o.inPlace ? o.in : o.out);
        map<const Genome *, GenomeData> data;
        map<const Genome *, string> likFiles;
        vector<const Genome *> pre = preorderAll(aln);
        for (auto it = pre.rbegin(); it != pre.rend(); ++it) {   /* children first */
            const Genome *A = *it;
            if (A->children().empty()) continue;
            vector<pair<const Genome *, Jc>> kids;
            for (const Genome *c : A->children()) kids.push_back({c, Jc(aln.manifest().find(c->name())->branchLength)});
            string file = scratch.file(genomeFileStem(A->name()) + ".hlk");
            map<string, string> meta{{"kind", "subtree"}, {"model", "JC69"}, {"tree", aln.manifest().newick},
                                     {"program", "hal likelihoods"}};
            vector<RefInterval> wins = windows(A);
            auto compute = [&](size_t k) {
                const RefInterval &w = wins[k];
                Worker wk(aln, data);
                int64_t n = w.end - w.start;
                vector<double> L((size_t)(4 * n), 1.0);
                vector<LiftPiece> in{LiftPiece{w.start, w.seq, w.start, n, false}}, out;
                for (auto &kc : kids) {
                    const Genome *c = kc.first;
                    const Jc &P = kc.second;
                    liftDownBranch(c, true, in, out, wk.cursors);
                    wk.forPieces(c, out, [&](const LiftPiece &pc, const double *v) {
                        for (int64_t i = 0; i < pc.length; i++) {
                            const double *cv = v + 4 * (pc.reversed ? pc.length - 1 - i : i);
                            double x[4];
                            for (int b = 0; b < 4; b++) x[b] = pc.reversed ? cv[b ^ 3] : cv[b];
                            double s = x[0] + x[1] + x[2] + x[3];
                            double *l = &L[(size_t)(4 * (pc.srcStart + i - w.start))];
                            for (int a = 0; a < 4; a++) l[a] *= P.other * s + (P.same - P.other) * x[a];
                            double t = l[0] + l[1] + l[2] + l[3];
                            if (t < 1e-200 && t > 0)
                                for (int a = 0; a < 4; a++) l[a] /= t;
                        }
                    });
                }
                WindowOut r;
                r.codes.resize((size_t)n);
                const string &b = wk.fetch(A, w.seq, w.start, w.end);
                for (int64_t i = 0; i < n; i++) {
                    double *l = &L[(size_t)(4 * i)];
                    normalize(l);
                    r.codes[(size_t)i] = encodeLikelihood(l, likBase(b[(size_t)i]), o.q);
                }
                return r;
            };
            runGenome(A, wins, file, meta, o, compute, nullptr, "");
            likFiles[A] = file;
            data[A].lik = LikFile::open(lmap_io_open_file(file.c_str()), true, file);
            cerr << "hal likelihoods: " << A->name() << endl;
        }
        attach(o, manifestBefore, data, likFiles, "likelihoods");
        return 0;
    });
}

static const char *POSTERIORS_USAGE =
    "usage: hal posteriors [options] <in.hal> <out.hal>\n"
    "       hal posteriors --inPlace [options] <in.hal>\n"
    "The probability of each base of every ancestral genome given all the data: a top-down pass\n"
    "over the ancestors' subtree likelihoods (from hal likelihoods, for the same tree), under\n"
    "the same model.  Each ancestral base is set to the most probable one, keeping its case\n"
    "(unless --keepBases), and the posteriors replace the likelihoods (kind \"posterior\").\n"
    "options:\n"
    "  --keepBases          leave the ancestral DNA as it is\n"
    "  --threads N          windows of a genome to do at once (default 1)\n"
    "  --stepsPerDecade Q   precision of the stored values: log10 ratios to 1/Q, 1 to 25 (default\n"
    "                       10); a file spans 254/Q decades, smaller probabilities stored as that\n"
    "  --inPlace            change in.hal itself (one crash-safe transaction for a single file)\n"
    "  --directory          write out.hal as a directory\n"
    "  -h, --help\n";

int main_posteriors(int argc, char **argv) {
    return runCommand("posteriors", [&]() {
        Options o;
        parse(argc, argv, POSTERIORS_USAGE, o, true);
        Alignment aln(o.in);
        string manifestBefore = aln.store().read(MANIFEST_NAME);
        preflight(o, aln);
        /* the subtree likelihoods every ancestor must have, for this tree, read with its
         * old DNA */
        map<const Genome *, GenomeData> old;
        for (const Genome *g : preorderAll(aln)) {
            if (g->children().empty()) continue;
            const MemberRef *mr = aln.manifest().find(g->name())->member("lik");
            if (!mr) throw runtime_error(g->name() + " has no likelihoods: run hal likelihoods first");
            unique_ptr<LikFile> f = LikFile::open(aln.store().openMember(mr->path), true, mr->path);
            if (f->meta("kind") != "subtree")
                throw runtime_error(g->name() + "'s likelihoods are \"" + f->meta("kind") + "\", not subtree likelihoods: run hal likelihoods first");
            if (f->meta("tree") != aln.manifest().newick)
                throw runtime_error(g->name() + "'s likelihoods were computed for another tree: run hal likelihoods again");
            if (!sameSequences(*f, g)) throw runtime_error(g->name() + "'s likelihoods are not for its sequences");
            old[g].lik = move(f);
        }
        Scratch scratch(o.inPlace ? o.in : o.out);
        map<const Genome *, GenomeData> post;       /* each done ancestor's posteriors (and new DNA) */
        map<const Genome *, string> likFiles;
        for (const Genome *c : preorderAll(aln)) {   /* parents first */
            if (c->children().empty()) continue;
            const Genome *p = c->parent();
            unique_ptr<Jc> P(p ? new Jc(aln.manifest().find(c->name())->branchLength) : nullptr);
            string file = scratch.file(genomeFileStem(c->name()) + ".hlk");
            string dnaFile = o.keepBases ? "" : scratch.file(genomeFileStem(c->name()) + ".2bit");
            map<string, string> meta{{"kind", "posterior"}, {"model", "JC69"}, {"tree", aln.manifest().newick},
                                     {"program", "hal posteriors"}};
            vector<string> names;
            for (uint32_t i = 0; i < c->numSequences(); i++) names.push_back(c->sequenceName(i));
            vector<const char *> cnames;
            for (const string &nm : names) cnames.push_back(nm.c_str());
            tb_writer *tw = nullptr;
            if (!o.keepBases) {
                tw = tb_writer_open(dnaFile.c_str(), cnames.data(), (uint32_t)cnames.size(), 0);
                if (!tw) throw runtime_error("cannot write " + dnaFile + ": " + strerror(errno));
            }
            unique_ptr<tb_writer, void (*)(tb_writer *)> twGuard(tw, [](tb_writer *w) {
                if (w) tb_writer_close(w, nullptr, 0);
            });
            vector<RefInterval> wins = windows(c);
            auto compute = [&](size_t k) {
                const RefInterval &w = wins[k];
                int64_t n = w.end - w.start;
                Worker mine(aln, old);    /* c's stored likelihoods, c's old DNA */
                Worker parent(aln, post); /* p's posteriors, p's new DNA */
                vector<double> L, U((size_t)(4 * n), 0.25);
                mine.vectors(c, w.seq, w.start, w.end, L);
                if (p) {
                    vector<LiftPiece> in{LiftPiece{w.start, w.seq, w.start, n, false}}, out;
                    liftUpBranch(c, in, out, nullptr, mine.cursors);
                    parent.forPieces(p, out, [&](const LiftPiece &pc, const double *pv) {
                        for (int64_t i = 0; i < pc.length; i++) {
                            int64_t ci = pc.srcStart + i - w.start;   /* the child base */
                            const double *pp = pv + 4 * (pc.reversed ? pc.length - 1 - i : i);
                            const double *lc = &L[(size_t)(4 * ci)];
                            /* both in the parent base's orientation */
                            double x[4], D[4], sd = 0;
                            for (int b = 0; b < 4; b++) x[b] = pc.reversed ? lc[b ^ 3] : lc[b];
                            double s = x[0] + x[1] + x[2] + x[3];
                            for (int a = 0; a < 4; a++) {
                                double m = P->other * s + (P->same - P->other) * x[a];
                                D[a] = m > 1e-300 ? pp[a] / m : 0;
                                sd += D[a];
                            }
                            double u[4];
                            for (int b = 0; b < 4; b++) {
                                u[b] = 0;
                                for (int a = 0; a < 4; a++) u[b] += (sd > 0 ? D[a] / sd : 0.25) * (a == b ? P->same : P->other);
                            }
                            double *uc = &U[(size_t)(4 * ci)];
                            for (int b = 0; b < 4; b++) uc[b] = pc.reversed ? u[b ^ 3] : u[b];
                        }
                    });
                }
                WindowOut r;
                r.bases = mine.fetch(c, w.seq, w.start, w.end);
                r.codes.resize((size_t)n);
                for (int64_t i = 0; i < n; i++) {
                    double v[4];
                    for (int b = 0; b < 4; b++) v[b] = L[(size_t)(4 * i + b)] * U[(size_t)(4 * i + b)];
                    normalize(v);
                    char &base = r.bases[(size_t)i];
                    if (!o.keepBases) {
                        double mx = max(max(v[0], v[1]), max(v[2], v[3]));
                        if (mx > 0.25 + 1e-12) {
                            int best = likBase(base) < 4 && v[likBase(base)] == mx ? likBase(base) : 0;
                            while (v[best] != mx) best++;
                            char nb = "ACGT"[best];
                            base = islower((unsigned char)base) ? (char)tolower(nb) : nb;
                        }
                    }
                    r.codes[(size_t)i] = encodeLikelihood(v, likBase(base), o.q);
                }
                return r;
            };
            runGenome(c, wins, file, meta, o, compute, tw, dnaFile);
            if (tw) {
                char err[256] = "";
                twGuard.release();
                if (tb_writer_close(tw, err, sizeof err) != 0) throw runtime_error(dnaFile + ": " + err);
                post[c].dnaFile = dnaFile;
            }
            likFiles[c] = file;
            post[c].lik = LikFile::open(lmap_io_open_file(file.c_str()), true, file);
            cerr << "hal posteriors: " << c->name() << endl;
        }
        attach(o, manifestBefore, post, likFiles, "posteriors");
        return 0;
    });
}
}
