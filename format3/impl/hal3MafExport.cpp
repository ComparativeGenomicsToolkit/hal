/*
 * HAL format 3: MAF and TAF export, row by row.  See hal3MafExport.h.
 */
#include "hal3MafExport.h"
#include "hal3Manifest.h"
#include "hal3Universal.h"
#include "twobit64.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <zlib.h>

using namespace std;

namespace hal3 {

/* ---------------------------------------------------------------- bytes out */

namespace {

class ByteSink {
  public:
    virtual ~ByteSink() {}
    virtual void write(const char *p, size_t n) = 0;
};

/* The output: a file, or stdout. */
class FileSink : public ByteSink {
  public:
    explicit FileSink(const string &path) : _own(path != "stdout") {
        _f = _own ? fopen(path.c_str(), "wb") : stdout;
        if (!_f) throw runtime_error("cannot write " + path + ": " + strerror(errno));
    }
    ~FileSink() override {
        if (_f && _own) fclose(_f);
    }
    void write(const char *p, size_t n) override {
        if (n && fwrite(p, 1, n, _f) != n) throw runtime_error("write failed");
    }
    void close() {
        FILE *f = _f;
        _f = nullptr;
        if (f && (_own ? fclose(f) != 0 : fflush(f) != 0)) throw runtime_error("write failed");
    }

  private:
    FILE *_f;
    bool _own;
};

/* One chunk's output, held until the chunks before it are written: in memory, or past
 * SPILL bytes in an unlinked temporary file -- in the first of tmpDirs where one can be
 * made. */
class ChunkBuffer : public ByteSink {
  public:
    explicit ChunkBuffer(const vector<string> &tmpDirs) : _tmpDirs(tmpDirs) {}
    ~ChunkBuffer() override {
        if (_tmp) fclose(_tmp);
    }
    ChunkBuffer(const ChunkBuffer &) = delete;
    ChunkBuffer &operator=(const ChunkBuffer &) = delete;
    void write(const char *p, size_t n) override {
        if (!_tmp && _mem.size() + n > SPILL) spill();
        if (!_tmp) _mem.append(p, n);
        else if (fwrite(p, 1, n, _tmp) != n) throw runtime_error("cannot write a temporary file in " + _tmpDir);
    }
    void copyTo(ByteSink &out) {
        if (!_tmp) {
            out.write(_mem.data(), _mem.size());
            return;
        }
        if (fflush(_tmp) != 0 || fseek(_tmp, 0, SEEK_SET) != 0) throw runtime_error("cannot read a temporary file in " + _tmpDir);
        vector<char> buf((size_t)1 << 20);
        size_t k;
        while ((k = fread(buf.data(), 1, buf.size(), _tmp)) > 0) out.write(buf.data(), k);
        if (ferror(_tmp)) throw runtime_error("cannot read a temporary file in " + _tmpDir);
    }

  private:
    static const size_t SPILL = (size_t)32 << 20;
    const vector<string> &_tmpDirs;
    string _tmpDir, _mem;
    FILE *_tmp = nullptr;
    void spill() {
        int fd = -1;
        string why;
        for (const string &d : _tmpDirs) {
            string t = d + "/.hal-maf-XXXXXX";
            fd = mkstemp(&t[0]);
            if (fd >= 0) {
                unlink(t.c_str());
                _tmpDir = d;
                break;
            }
            why += (why.empty() ? "" : ", ") + d + ": " + strerror(errno);
        }
        if (fd < 0) throw runtime_error("cannot create a temporary file (" + why + ")");
        _tmp = fdopen(fd, "w+b");
        if (!_tmp) {
            ::close(fd);
            throw runtime_error("cannot create a temporary file in " + _tmpDir);
        }
        if (!_mem.empty() && fwrite(_mem.data(), 1, _mem.size(), _tmp) != _mem.size())
            throw runtime_error("cannot write a temporary file in " + _tmpDir);
        string().swap(_mem);
    }
};

/* Text, written plain or as BGZF: gzip members of at most 64 KiB, each with a BC extra
 * field giving its size, as bgzip writes.  A block is cut every 65280 bytes of text and
 * at finish(), which adds no end-of-file marker (see writeBgzfEof). */
class TextOut {
  public:
    TextOut(ByteSink &out, bool bgzf) : _out(out), _bgzf(bgzf) {}
    ~TextOut() {
        if (_zInit) deflateEnd(&_z);
    }
    TextOut(const TextOut &) = delete;
    TextOut &operator=(const TextOut &) = delete;
    void put(const char *p, size_t n) {
        if (!_bgzf) {
            _out.write(p, n);
            return;
        }
        while (n) {
            size_t k = min(n, BLOCK - _buf.size());
            _buf.append(p, k);
            p += k;
            n -= k;
            if (_buf.size() == BLOCK) block();
        }
    }
    void put(const string &s) { put(s.data(), s.size()); }
    void finish() {
        if (_bgzf && !_buf.empty()) block();
    }

  private:
    static const size_t BLOCK = 65280;
    ByteSink &_out;
    bool _bgzf;
    string _buf;
    z_stream _z;
    bool _zInit = false;
    void block() {
        if (!_zInit) {
            memset(&_z, 0, sizeof _z);
            if (deflateInit2(&_z, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
                throw runtime_error("deflateInit2 failed");
            _zInit = true;
        } else if (deflateReset(&_z) != Z_OK) {
            throw runtime_error("deflateReset failed");
        }
        unsigned char out[65536];
        size_t n = _buf.size();
        _z.next_in = (Bytef *)_buf.data();
        _z.avail_in = (uInt)n;
        _z.next_out = out + 18;
        _z.avail_out = sizeof out - 18 - 8;
        int rc = deflate(&_z, Z_FINISH);
        size_t clen = _z.total_out;
        if (rc != Z_STREAM_END) throw runtime_error("BGZF block did not fit");
        size_t total = 18 + clen + 8;
        unsigned char hdr[18] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0,
                                 (unsigned char)((total - 1) & 0xff), (unsigned char)((total - 1) >> 8)};
        memcpy(out, hdr, 18);
        uLong crc = crc32(crc32(0L, Z_NULL, 0), (const Bytef *)_buf.data(), (uInt)n);
        unsigned char *t = out + 18 + clen;
        for (int k = 0; k < 4; k++) t[k] = (unsigned char)(crc >> (8 * k));
        for (int k = 0; k < 4; k++) t[4 + k] = (unsigned char)(n >> (8 * k));
        _out.write((const char *)out, total);
        _buf.clear();
    }
};

void writeBgzfEof(ByteSink &out) {
    static const unsigned char eof[28] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43,
                                          2,    0,    0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    out.write((const char *)eof, sizeof eof);
}

/* ---------------------------------------------------------------- blocks */

/* A finished block's row.  id is the row's chain (nullptr for row 0, the reference or
 * the universal top-level range): the row continues the previous block's row with the
 * same id. */
struct BlockRow {
    const void *id;
    const Genome *genome;
    uint32_t seq;
    bool reversed;
    int64_t start, size;   /* oriented (MAF) coordinates */
    string text;
};

class BlockWriter {
  public:
    virtual ~BlockWriter() {}
    /* rows[0 .. n) with rows[0] row 0 (see BlockRow); every text ncol long */
    virtual void block(const vector<BlockRow> &rows, size_t n, int64_t ncol) = 0;
};

string rowName(const Genome *g, uint32_t seq, bool onlySequenceNames) {
    return onlySequenceNames ? string(g->sequenceName(seq)) : g->name() + "." + g->sequenceName(seq);
}

class MafBlockWriter : public BlockWriter {
  public:
    MafBlockWriter(TextOut &out, bool onlySequenceNames) : _out(out), _names(onlySequenceNames) {}
    void block(const vector<BlockRow> &rows, size_t n, int64_t) override {
        _text.assign("a\n");
        for (size_t i = 0; i < n; i++) {
            const BlockRow &r = rows[i];
            _text += "s\t";
            _text += rowName(r.genome, r.seq, _names);
            _text += "\t" + to_string(r.start) + "\t" + to_string(r.size) + "\t" + (r.reversed ? "-" : "+") + "\t" +
                     to_string(r.genome->sequenceLength(r.seq)) + "\t";
            _text += r.text;
            _text += "\n";
        }
        _text += "\n";
        _out.put(_text);
    }

  private:
    TextOut &_out;
    bool _names;
    string _text;
};

/* TAF: a row's coordinates, as an i or s operation gives them. */
struct TafAnchor {
    string name;
    int64_t start;
    bool reversed;
    int64_t seqLength;
};

void appendAnchor(string &line, char op, size_t row, const TafAnchor &a) {
    line += ' ';
    line += op;
    line += ' ';
    line += to_string(row);
    line += ' ';
    line += a.name;
    line += ' ';
    line += to_string(a.start);
    line += a.reversed ? " - " : " + ";
    line += to_string(a.seqLength);
}

void appendCount(string &line, size_t k) {
    char buf[24];
    int i = 24;
    do {
        buf[--i] = (char)('0' + k % 10);
        k /= 10;
    } while (k);
    line.append(buf + i, (size_t)(24 - i));
}

/* One column's bases, col[0 .. n), as taffy writes them: each base, or with runLength
 * runs of equal bases as "base count ". */
void appendColumn(const char *col, size_t n, bool runLength, string &line) {
    if (!runLength) {
        line.append(col, n);
        return;
    }
    for (size_t i = 0; i < n;) {
        size_t j = i + 1;
        while (j < n && col[j] == col[i]) j++;
        line += col[i];
        line += ' ';
        appendCount(line, j - i);
        line += ' ';
        i = j;
    }
}

/* What one chunk produced. */
struct ChunkResult {
    unique_ptr<ChunkBuffer> body;
    MafStats stats;
    /* universal: the chunk's axis, and the moments of what each genome's rows hold */
    vector<Segment> axis;
    map<const Genome *, Moments> moments;
    /* TAF: the chunk's first line depends on how many rows the chunk before it ended
     * with, so it is left out of body: its first column and its rows' coordinates. */
    string firstColumn;
    vector<TafAnchor> firstRows;
    size_t lastRows = 0;
};

/* The TAF of one chunk, but for its first line (see ChunkResult).  Rows continue from
 * block to block by id; failing that, a row that ended continues into an unmatched row
 * of the same sequence and strand that starts at or after its end ("g" skips the bases
 * between).  Links that would reorder rows are dropped (a longest increasing run is
 * kept), as TAF cannot move a row.  Coordinates are repeated as taffy view repeats
 * them: for a row more than REPEAT bases since its last, and for every row whenever
 * row 0's are given.  columnCadence (universal, where row 0 may have no bases for
 * millions of columns): row 0's too once REPEAT columns have passed since every row's
 * were last given, so taffy's index finds coordinates for all rows that often. */
class TafBlockWriter : public BlockWriter {
  public:
    TafBlockWriter(TextOut &out, bool runLength, bool onlySequenceNames, bool columnCadence, ChunkResult &result)
        : _out(out), _runLength(runLength), _names(onlySequenceNames), _cadence(columnCadence), _result(result) {}
    void block(const vector<BlockRow> &rows, size_t n, int64_t ncol) override {
        _next.resize(n);
        for (size_t j = 0; j < n; j++) {
            const BlockRow &r = rows[j];
            _next[j] = Row{r.id, r.genome, r.seq, r.reversed, r.start, r.size, 0};
        }
        _line.clear();
        /* the block's columns, TILE at a time, transposed so each column is contiguous */
        for (int64_t c0 = 0; c0 < ncol; c0 += TILE) {
            size_t w = (size_t)min<int64_t>(TILE, ncol - c0);
            _cols.resize(w * n);
            for (size_t j = 0; j < n; j++) {
                const char *t = rows[j].text.data() + c0;
                for (size_t c = 0; c < w; c++) _cols[c * n + j] = t[c];
            }
            for (size_t c = 0; c < w; c++) {
                const char *col = _cols.data() + c * n;
                if (c0 + (int64_t)c > 0) {
                    appendColumn(col, n, _runLength, _line);
                } else if (_first) {
                    _first = false;
                    appendColumn(col, n, _runLength, _result.firstColumn);
                    for (size_t j = 0; j < n; j++) _result.firstRows.push_back(anchor(_next[j]));
                    if (_onFirst) _onFirst();
                    continue;
                } else {
                    appendColumn(col, n, _runLength, _line);
                    operations();
                }
                _line += '\n';
            }
            if (_line.size() >= ((size_t)1 << 16)) {
                _out.put(_line);
                _line.clear();
            }
        }
        _out.put(_line);
        _prev.swap(_next);
        _colsSinceAll += ncol;
    }
    void finish() { _result.lastRows = _prev.size(); }
    /* called once the chunk's first line is known (its first column and rows), before
     * anything else of the chunk is written */
    void onFirstLine(const function<void()> &f) { _onFirst = f; }

  private:
    static const int64_t REPEAT = 10000;   /* taffy view's default */
    static const int64_t TILE = 256;
    struct Row {
        const void *id;
        const Genome *genome;
        uint32_t seq;
        bool reversed;
        int64_t start, size, since;   /* since: bases since its coordinates were last given */
    };
    TextOut &_out;
    bool _runLength, _names, _cadence;
    ChunkResult &_result;
    int64_t _colsSinceAll = 0;   /* columns since every row's coordinates were given */
    bool _first = true;
    function<void()> _onFirst;
    vector<Row> _prev, _next;
    string _line, _cols;

    TafAnchor anchor(const Row &r) const {
        return TafAnchor{rowName(r.genome, r.seq, _names), r.start, r.reversed, r.genome->sequenceLength(r.seq)};
    }
    /* _line += the coordinate operations taking _prev's rows to _next's */
    void operations() {
        size_t p = _prev.size(), n = _next.size();
        vector<int64_t> link(n, -1);   /* next row -> the prev row it continues */
        link[0] = 0;                   /* row 0: always, and always adjacent (one sequence per chunk) */
        unordered_map<const void *, size_t> byId;
        for (size_t i = 1; i < p; i++) byId[_prev[i].id] = i;
        vector<bool> used(p, false);
        used[0] = true;
        for (size_t j = 1; j < n; j++) {
            auto f = byId.find(_next[j].id);
            if (f != byId.end() && follows(_prev[f->second], _next[j])) {
                link[j] = (int64_t)f->second;
                used[f->second] = true;
            }
        }
        /* rows that ended, continued by new rows of the same sequence and strand */
        map<tuple<const Genome *, uint32_t, bool>, vector<size_t>> ended;
        for (size_t i = 1; i < p; i++)
            if (!used[i]) ended[make_tuple(_prev[i].genome, _prev[i].seq, _prev[i].reversed)].push_back(i);
        if (!ended.empty()) {
            map<tuple<const Genome *, uint32_t, bool>, size_t> from;   /* where each list's search resumes */
            for (size_t j = 1; j < n; j++) {
                if (link[j] >= 0) continue;
                auto key = make_tuple(_next[j].genome, _next[j].seq, _next[j].reversed);
                auto e = ended.find(key);
                if (e == ended.end()) continue;
                size_t &k = from[key];
                for (size_t q = k; q < e->second.size(); q++) {
                    size_t i = e->second[q];
                    if (!used[i] && follows(_prev[i], _next[j])) {
                        link[j] = (int64_t)i;
                        used[i] = true;
                        k = q + 1;
                        break;
                    }
                }
            }
        }
        keepIncreasing(link);
        fill(used.begin(), used.end(), false);
        for (size_t j = 0; j < n; j++)
            if (link[j] >= 0) used[(size_t)link[j]] = true;
        /* deletions first, each numbered among the rows still there */
        _line += " ;";
        size_t kept = 0;
        for (size_t i = 0; i < p; i++) {
            if (used[i]) {
                kept++;
                continue;
            }
            _line += " d ";
            _line += to_string(kept);
        }
        /* then insertions, coordinates and skips, in row order */
        bool everything = false;
        for (size_t j = 0; j < n; j++) {
            Row &r = _next[j];
            if (link[j] < 0) {
                appendAnchor(_line, 'i', j, anchor(r));
                r.since = 0;
                if (j == 0) everything = true;
                continue;
            }
            const Row &q = _prev[(size_t)link[j]];
            r.since = q.since + q.size;
            if (everything || r.since > REPEAT || (j == 0 && _cadence && _colsSinceAll >= REPEAT)) {
                appendAnchor(_line, 's', j, anchor(r));
                r.since = 0;
                if (j == 0) everything = true;
                continue;
            }
            int64_t gap = r.start - (q.start + q.size);
            if (gap > 0) {
                _line += " g ";
                _line += to_string(j);
                _line += ' ';
                _line += to_string(gap);
            }
        }
        if (everything) _colsSinceAll = 0;
    }
    static bool follows(const Row &q, const Row &r) {
        return q.genome == r.genome && q.seq == r.seq && q.reversed == r.reversed && r.start >= q.start + q.size;
    }
    /* keep the longest run of links increasing in both rows (the reference row's
     * first); unlink the rest */
    static void keepIncreasing(vector<int64_t> &link) {
        vector<size_t> js;
        for (size_t j = 1; j < link.size(); j++)
            if (link[j] >= 0) js.push_back(j);
        vector<size_t> tails;               /* tails[l]: index into js ending a run of length l+1 */
        vector<int64_t> before(js.size(), -1);
        for (size_t t = 0; t < js.size(); t++) {
            int64_t v = link[js[t]];
            size_t lo = 0, hi = tails.size();
            while (lo < hi) {
                size_t mid = (lo + hi) / 2;
                if (link[js[tails[mid]]] < v) lo = mid + 1;
                else hi = mid;
            }
            if (lo > 0) before[t] = (int64_t)tails[lo - 1];
            if (lo == tails.size()) tails.push_back(t);
            else tails[lo] = t;
        }
        vector<bool> keep(js.size(), false);
        for (int64_t t = tails.empty() ? -1 : (int64_t)tails.back(); t >= 0; t = before[(size_t)t]) keep[(size_t)t] = true;
        for (size_t t = 0; t < js.size(); t++)
            if (!keep[t]) link[js[t]] = -1;
    }
};

/* The line starting a chunk, after a chunk that ended with prevRows rows: its first
 * column, rows past the chunk's own deleted, then every row given coordinates -- in
 * place (s) where a row is there already, else inserted (i). */
string tafFirstLine(const ChunkResult &r, size_t prevRows) {
    string line = r.firstColumn + " ;";
    size_t n = r.firstRows.size();
    for (size_t k = n; k < prevRows; k++) line += " d " + to_string(n);
    for (size_t j = 0; j < n; j++) appendAnchor(line, j < prevRows ? 's' : 'i', j, r.firstRows[j]);
    return line + "\n";
}

/* Genomes in row order: by distance from the reference in the tree's preorder, ties
 * the earlier (cactus-hal2maf's genome list); any genome the tree does not name last. */
map<const Genome *, int> rowRanks(const Alignment &aln, const Genome *ref) {
    vector<const Genome *> pre;
    set<const Genome *> seen;
    for (const NewickNode &n : parseNewick(aln.manifest().newick)) {
        const Genome *g = aln.genome(n.name);
        if (g && seen.insert(g).second) pre.push_back(g);
    }
    int at = (int)(find(pre.begin(), pre.end(), ref) - pre.begin());
    vector<int> order((size_t)pre.size());
    for (size_t i = 0; i < order.size(); i++) order[i] = (int)i;
    stable_sort(order.begin(), order.end(), [&](int a, int b) { return abs(a - at) < abs(b - at); });
    map<const Genome *, int> rank;
    for (size_t i = 0; i < order.size(); i++) rank[pre[(size_t)order[i]]] = (int)i;
    for (const GenomeEntry &e : aln.manifest().genomes) {
        const Genome *g = aln.genome(e.name);
        if (g && !rank.count(g)) rank[g] = (int)rank.size();
    }
    return rank;
}

WalkOptions walkOptions(const MafOptions &o) {
    WalkOptions w;
    w.ancestors = o.ancestors;
    w.dupes = o.dupes;
    w.universal = o.universal;
    w.allAncestorRows = o.allAncestorRows;
    w.onlyOrthologs = o.onlyOrthologs;
    w.targets = o.targets;
    w.rootGenome = o.rootGenome;
    w.exclude = o.exclude;
    w.refSingleCopy = o.refSingleCopy;
    return w;
}

/* ---------------------------------------------------------------- one chunk */

/* A universal chunk: its segments (hal3Universal.h), and its top-level range's sequence,
 * whose bases from rangePos are the chunk's depth-0 segments. */
struct UniChunk {
    const UniversalIndex *index;
    const vector<Segment> *segs;
    const Genome *genome;
    uint32_t seq;
    int64_t rangePos, cols;
};

/* Exports chunks one at a time; one per thread, kept from chunk to chunk with its walk
 * (so its cursors) and its DNA handles. */
class ChunkExporter {
  public:
    ChunkExporter(const Alignment &alignment, const Genome *ref, const MafOptions &opts, const map<const Genome *, int> &ranks)
        : _aln(alignment), _ref(ref), _opts(opts), _ranks(ranks), _walk(alignment, ref, walkOptions(opts)), _dna(alignment) {}
    ChunkExporter(const ChunkExporter &) = delete;
    ChunkExporter &operator=(const ChunkExporter &) = delete;

    void run(const RefInterval &chunk, BlockWriter &out);
    void runUniversal(const UniChunk &chunk, BlockWriter &out);
    MafStats stats;                          /* of the last run */
    map<const Genome *, Moments> moments;    /* of the last universal run (universalWhole) */

  private:
    /* Row 0: the reference, at every column (column x = position x); or the universal
     * range's sequence, its bases from pos at the columns runs [u, u+len) in order. */
    struct RowZero {
        const Genome *genome;
        uint32_t seq;
        bool all;
        int64_t pos;
        vector<int64_t> runU, runLen, before;   /* before[i]: its bases in runs before i */
        void add(int64_t u, int64_t len) {
            before.push_back(before.empty() ? 0 : before.back() + runLen.back());
            runU.push_back(u);
            runLen.push_back(len);
        }
        /* its position at column x: the next base's */
        int64_t at(int64_t x) const {
            if (all) return x;
            size_t i = (size_t)(upper_bound(runU.begin(), runU.end(), x) - runU.begin());
            if (i == 0) return pos;
            i--;
            return pos + before[i] + min(runLen[i], x - runU[i]);
        }
    };
    struct Key {   /* a top ancestral base: where a base's homology is decided */
        const Genome *node;
        uint32_t seq;
        int64_t pos;
        bool operator<(const Key &o) const {
            if (node != o.node) return node->name() < o.node->name();
            if (seq != o.seq) return seq < o.seq;
            return pos < o.pos;
        }
        bool operator==(const Key &o) const { return node == o.node && seq == o.seq && pos == o.pos; }
    };
    /* A row: pieces of one genome, contiguous in it except for inserted bases between
     * them, in reference order.  Oriented coordinates (MAF's: position on the row's
     * strand) increase along the reference.  Link i joins pieces i and i+1; the bases
     * between them, [oend(i), ostart[i+1]), are inserted: they go in the insertion group
     * at reference boundary at[i] (between reference positions at[i]-1 and at[i]). */
    struct Chain {
        const Genome *genome;
        uint32_t seq;
        bool reversed;
        vector<LiftPiece> pieces;   /* in reference order */
        vector<int64_t> ostart;     /* oriented start of each piece */
        vector<int64_t> at;         /* per link: its boundary, or -1 if nothing is inserted */
        vector<vector<Key>> keys;   /* per link: each inserted base's top key */
        vector<vector<int>> cols;   /* per link: each inserted base's column in its group */
        int64_t start() const { return pieces.front().srcStart; }
        int64_t end() const { return pieces.back().srcStart + pieces.back().length; }
        int64_t oend(size_t i) const { return ostart[i] + pieces[i].length; }
    };
    const Alignment &_aln;
    const Genome *_ref;
    const MafOptions &_opts;
    const map<const Genome *, int> &_ranks;
    TreeWalk _walk;
    DnaReader _dna;   /* recycled at each chunk: what its handles hold stays about one chunk's */

    int rank(const Genome *g) const {
        auto r = _ranks.find(g);
        return r == _ranks.end() ? INT_MAX : r->second;
    }
    /* output order of rows (after the reference row) */
    bool rowBefore(const Chain &a, int64_t ao, const Chain &b, int64_t bo) const {
        if (a.genome != b.genome) return rank(a.genome) < rank(b.genome);
        if (a.seq != b.seq) {
            int c = strcmp(a.genome->sequenceName(a.seq), b.genome->sequenceName(b.seq));
            if (c) return c < 0;
        }
        if (a.reversed != b.reversed) return !a.reversed;
        return ao < bo;
    }
    bool insertionKeys(const Genome *g, uint32_t seq, bool reversed, int64_t o1, int64_t o2, vector<Key> &keys);
    void startRun();
    /* blocks of columns [cs, ce) from the lifted pieces (sources: columns); cut
     * candidates are chain ends and extra */
    void layout(PieceMap &lifted, int64_t cs, int64_t ce, const RowZero &z, const vector<int64_t> &extra, BlockWriter &writer);
    void fetch(const Genome *g, uint32_t seq, int64_t start, int64_t end, bool reversed, string &out);
};

/* ---------------------------------------------------------------- insertions */

/* The bases [o1, o2) of g's sequence seq, in oriented coordinates (on the reverse strand
 * when reversed), as candidate inserted bases.  keys gets each base's top ancestral base,
 * in oriented order.  Returns false if any of them is in a column holding a reference
 * base -- such a base belongs in that column, not in an insertion.
 *
 * Default and --noDupes: a base's column is the down-closure of its top base (walking up
 * with the mode's copy rule), so it holds a reference base exactly when walking down from
 * the top towards the reference reaches one.  --onlyOrthologs: when liftover from the
 * base to the reference (through their common ancestor) reaches one. */
bool ChunkExporter::insertionKeys(const Genome *g, uint32_t seq, bool reversed, int64_t o1, int64_t o2, vector<Key> &keys) {
    int64_t len = o2 - o1, seqLen = g->sequenceLength(seq);
    int64_t lo = reversed ? seqLen - o2 : o1;   /* forward coordinates */
    keys.assign((size_t)len, Key{nullptr, 0, 0});
    if (_opts.onlyOrthologs) {
        vector<LiftPiece> hits;
        liftInterval(g, seq, lo, lo + len, _ref, _opts.dupes, hits, &_walk.cursors());
        if (!hits.empty()) return false;
    }
    /* up to the tops */
    vector<LiftPiece> cur{LiftPiece{lo, seq, lo, len, false}}, next, unmapped;
    vector<pair<const Genome *, vector<LiftPiece>>> tops;
    for (const Genome *a = g; !cur.empty(); a = a->parent()) {
        if (a == _walk.top()) {
            tops.emplace_back(a, cur);
            break;
        }
        liftUpBranch(a, cur, next, &unmapped, _walk.cursors(), !_opts.dupes);
        if (!unmapped.empty()) tops.emplace_back(a, unmapped);
        cur.swap(next);
    }
    for (auto &t : tops) {
        for (const LiftPiece &p : t.second) {
            for (int64_t i = 0; i < p.length; i++) {
                int64_t fwd = p.srcStart + i;                    /* the base in g, forward */
                int64_t at = p.reversed ? p.tgtStart + p.length - 1 - i : p.tgtStart + i;
                int64_t k = reversed ? (lo + len - 1 - fwd) : (fwd - lo);   /* oriented index */
                keys[(size_t)k] = Key{t.first, p.tgtSeq, at};
            }
        }
        if (_opts.onlyOrthologs) continue;
        /* down from the top towards the reference, if the top is on its lineage */
        const vector<const Genome *> &path = _walk.refPath();
        auto on = find(path.begin(), path.end(), t.first);
        if (on == path.end()) continue;
        vector<LiftPiece> down = t.second;
        for (auto c = on; c != path.begin() && !down.empty();) {
            --c;
            liftDownBranch(*c, _opts.dupes, down, next, _walk.cursors());
            down.swap(next);
        }
        if (!down.empty()) return false;
    }
    /* bases with no key cannot happen (every base tops out somewhere) */
    for (const Key &k : keys)
        if (!k.node) throw runtime_error("internal error: an inserted base with no top");
    return true;
}

/* ---------------------------------------------------------------- rows and blocks */

void ChunkExporter::fetch(const Genome *g, uint32_t seq, int64_t start, int64_t end, bool reversed, string &out) {
    if (_opts.noDna) {
        out.assign((size_t)(end - start), 'N');
        return;
    }
    const GenomeEntry *e = _aln.manifest().find(g->name());
    if (!e || e->dnaPath.empty()) throw runtime_error("genome " + g->name() + " has no DNA in this alignment (--noDna writes N)");
    int64_t len = g->sequenceLength(seq);
    _dna.fetch(g, seq, reversed ? len - end : start, reversed ? len - start : end, out);   /* forward coordinates */
    if (reversed) reverseComplement(out);
}

void ChunkExporter::startRun() {
    stats = MafStats();
    _dna.recycle();
}

void ChunkExporter::run(const RefInterval &chunk, BlockWriter &writer) {
    startRun();
    PieceMap lifted;
    _walk.lift(chunk.seq, {{chunk.start, chunk.end}}, lifted);
    RowZero z{_ref, chunk.seq, true, 0, {}, {}, {}};
    layout(lifted, chunk.start, chunk.end, z, {}, writer);
}

void ChunkExporter::runUniversal(const UniChunk &c, BlockWriter &writer) {
    startRun();
    moments.clear();
    const UniversalIndex &idx = *c.index;
    PieceMap own, quiet, lifted;
    RowZero z{c.genome, c.seq, false, c.rangePos, {}, {}, {}};
    vector<int64_t> extra;
    for (const Segment &sg : *c.segs) {
        (sg.depth ? own : quiet)[idx.genome(sg.g)].push_back(LiftPiece{sg.u, sg.seq, sg.pos, sg.len, sg.rev});
        if (sg.depth == 0) z.add(sg.u, sg.len);
        if (sg.depth > 0 || c.genome != _ref) stats.placedColumns += (uint64_t)sg.len;   /* an extra range is placed */
        if (sg.u > 0) extra.push_back(sg.u);
    }
    const int64_t rowZeroBases = z.at(c.cols) - c.rangePos;
    if (c.genome == _ref) stats.backboneColumns = (uint64_t)rowZeroBases;
    _walk.down(own, quiet, lifted);
    if (_opts.universalWhole) {
        for (auto &gp : lifted) {
            Moments &m = moments[gp.first];
            for (const LiftPiece &p : gp.second) {
                int64_t a = gp.first->sequenceStart(p.tgtSeq) + p.tgtStart;
                m.add(a, a + p.length);
            }
        }
        int64_t a = c.genome->sequenceStart(c.seq) + c.rangePos;
        moments[c.genome].add(a, a + rowZeroBases);
    }
    layout(lifted, 0, c.cols, z, extra, writer);
}

void ChunkExporter::layout(PieceMap &lifted, const int64_t cs, const int64_t ce, const RowZero &z, const vector<int64_t> &extra,
                           BlockWriter &writer) {
    const int64_t n = _opts.universal ? 0 : max<int64_t>(0, _opts.joinGap);

    /* ---- rows: link each genome's pieces where one follows another in the genome --
     * directly, or across at most n bases that are in no reference column -- and later
     * on the reference.  Links are taken greedily: fewest inserted bases, then fewest
     * skipped reference bases. */
    vector<Chain> chains;
    for (auto &gp : lifted) {
        const Genome *g = gp.first;
        vector<LiftPiece> &ps = gp.second;
        size_t np = ps.size();
        vector<int64_t> os(np);
        for (size_t i = 0; i < np; i++) {
            const LiftPiece &p = ps[i];
            os[i] = p.reversed ? g->sequenceLength(p.tgtSeq) - (p.tgtStart + p.length) : p.tgtStart;
        }
        map<pair<uint32_t, bool>, vector<pair<int64_t, size_t>>> byStart;   /* oriented start -> piece */
        for (size_t i = 0; i < np; i++) byStart[make_pair(ps[i].tgtSeq, ps[i].reversed)].push_back({os[i], i});
        for (auto &b : byStart) sort(b.second.begin(), b.second.end());
        vector<tuple<int64_t, int64_t, int64_t, size_t, size_t>> edges;   /* dTgt, dRef, p.src, p, q */
        for (size_t i = 0; i < np; i++) {
            const vector<pair<int64_t, size_t>> &v = byStart[make_pair(ps[i].tgtSeq, ps[i].reversed)];
            int64_t oe = os[i] + ps[i].length, pend = ps[i].srcStart + ps[i].length;
            for (auto it = lower_bound(v.begin(), v.end(), make_pair(oe, (size_t)0)); it != v.end() && it->first <= oe + n; ++it)
                if (ps[it->second].srcStart >= pend) edges.emplace_back(it->first - oe, ps[it->second].srcStart - pend, ps[i].srcStart, i, it->second);
        }
        sort(edges.begin(), edges.end());
        vector<int64_t> succ(np, -1), pred(np, -1);
        vector<bool> blocked(np, false);   /* p's nearest gap holds a reference-aligned base */
        map<size_t, vector<Key>> linkKeys;  /* by p */
        vector<Key> keys;
        for (auto &e : edges) {
            int64_t dTgt = get<0>(e);
            size_t p = get<3>(e), q = get<4>(e);
            if (succ[p] >= 0 || pred[q] >= 0 || blocked[p]) continue;
            if (dTgt > 0) {
                /* a wider gap for p contains this one, so one refusal settles p */
                if (!insertionKeys(g, ps[p].tgtSeq, ps[p].reversed, os[p] + ps[p].length, os[q], keys)) {
                    blocked[p] = true;
                    stats.joinsRejected++;
                    continue;
                }
                linkKeys[p] = keys;
            }
            succ[p] = (int64_t)q;
            pred[q] = (int64_t)p;
        }
        for (size_t i = 0; i < np; i++) {
            if (pred[i] >= 0) continue;
            Chain c{g, ps[i].tgtSeq, ps[i].reversed, {}, {}, {}, {}, {}};
            for (int64_t k = (int64_t)i; k >= 0; k = succ[k]) {
                if (!c.pieces.empty()) {
                    auto lk = linkKeys.find((size_t)pred[k]);
                    c.at.push_back(lk == linkKeys.end() ? -1 : 0);
                    c.keys.push_back(lk == linkKeys.end() ? vector<Key>() : lk->second);
                    c.cols.emplace_back();
                }
                c.pieces.push_back(ps[k]);
                c.ostart.push_back(os[k]);
            }
            chains.push_back(move(c));
        }
    }
    lifted.clear();
    sort(chains.begin(), chains.end(), [](const Chain &a, const Chain &b) {
        if (a.start() != b.start()) return a.start() < b.start();
        int c = a.genome->name().compare(b.genome->name());
        if (c) return c < 0;
        return make_tuple(a.seq, a.ostart[0], a.reversed) < make_tuple(b.seq, b.ostart[0], b.reversed);
    });

    /* ---- insertion groups.  A link that skips no reference base goes at the boundary
     * where its pieces meet; one that also skips reference bases (a deletion as well)
     * goes at whichever end shares more top keys with what is placed there already
     * (ties: the left), with its deletion gaps after it. */
    struct Member {
        size_t chain, link;
    };
    struct Group {
        vector<Member> members;
        set<Key> keys;
        int width = 0;
        int64_t bases = 0;
    };
    map<int64_t, Group> groups;
    auto place = [&](size_t ci, size_t li, int64_t b) {
        Chain &c = chains[ci];
        c.at[li] = b;
        Group &gr = groups[b];
        gr.members.push_back(Member{ci, li});
        gr.keys.insert(c.keys[li].begin(), c.keys[li].end());
        gr.bases += (int64_t)c.keys[li].size();
    };
    for (int pass = 0; pass < 2; pass++) {
        for (size_t ci = 0; ci < chains.size(); ci++) {
            Chain &c = chains[ci];
            for (size_t li = 0; li < c.at.size(); li++) {
                if (c.at[li] < 0) continue;
                int64_t left = c.pieces[li].srcStart + c.pieces[li].length, right = c.pieces[li + 1].srcStart;
                if ((pass == 0) != (left == right)) continue;
                if (left == right) {
                    place(ci, li, left);
                    continue;
                }
                auto shared = [&](int64_t b) {
                    auto g = groups.find(b);
                    if (g == groups.end()) return 0;
                    int k = 0;
                    for (const Key &key : c.keys[li]) k += (int)g->second.keys.count(key);
                    return k;
                };
                place(ci, li, shared(right) > shared(left) ? right : left);
            }
        }
    }
    /* Layout: members in output order; each base joins the first column after the row's
     * previous base that holds its key, or gets a new column right after it. */
    for (auto &gb : groups) {
        Group &gr = gb.second;
        stable_sort(gr.members.begin(), gr.members.end(), [&](const Member &x, const Member &y) {
            const Chain &a = chains[x.chain], &b = chains[y.chain];
            return rowBefore(a, a.oend(x.link), b, b.oend(y.link));
        });
        vector<Key> colKey;
        vector<int> order, pos;
        map<Key, vector<int>> withKey;
        for (const Member &m : gr.members) {
            Chain &c = chains[m.chain];
            vector<int> &ids = c.cols[m.link];
            ids.clear();
            int prev = -1;
            for (const Key &k : c.keys[m.link]) {
                int found = -1;
                auto w = withKey.find(k);
                if (w != withKey.end())
                    for (int id : w->second)
                        if (pos[id] > prev && (found < 0 || pos[id] < pos[found])) found = id;
                if (found < 0) {
                    found = (int)colKey.size();
                    colKey.push_back(k);
                    order.insert(order.begin() + (prev + 1), found);
                    pos.resize(colKey.size());
                    for (size_t i = 0; i < order.size(); i++) pos[order[i]] = (int)i;
                    withKey[k].push_back(found);
                }
                ids.push_back(found);
                prev = pos[found];
            }
        }
        for (const Member &m : gr.members)
            for (int &id : chains[m.chain].cols[m.link]) id = pos[id];
        gr.width = (int)order.size();
    }

    /* ---- counts for the block rule: chain starts and ends, covered bases, groups */
    vector<int64_t> starts, ends;
    vector<pair<int64_t, int>> cov;
    for (const Chain &c : chains) {
        starts.push_back(c.start());
        ends.push_back(c.end());
        for (const LiftPiece &p : c.pieces) {
            cov.push_back({p.srcStart, 1});
            cov.push_back({p.srcStart + p.length, -1});
        }
    }
    sort(starts.begin(), starts.end());
    sort(ends.begin(), ends.end());
    sort(cov.begin(), cov.end());
    vector<int64_t> covPos{cs}, covCum{0};   /* aligned reference-column bases before covPos[i] */
    vector<int> covDepth{0};
    {
        int depth = 0;
        for (auto &e : cov) {
            if (e.first != covPos.back()) {
                covCum.push_back(covCum.back() + (int64_t)depth * (e.first - covPos.back()));
                covPos.push_back(e.first);
                covDepth.push_back(depth);
            }
            depth += e.second;
            covDepth.back() = depth;
        }
    }
    vector<pair<int64_t, int>>().swap(cov);
    auto aligned = [&](int64_t x) {   /* aligned bases in reference columns [cs, x) */
        size_t i = upper_bound(covPos.begin(), covPos.end(), x) - covPos.begin() - 1;
        return covCum[i] + (int64_t)covDepth[i] * (x - covPos[i]);
    };
    vector<int64_t> groupAt{cs}, groupCols{0}, groupBases{0};   /* prefix sums over boundaries */
    for (auto &gb : groups) {
        groupAt.push_back(gb.first);
        groupCols.push_back(groupCols.back() + gb.second.width);
        groupBases.push_back(groupBases.back() + gb.second.bases);
    }
    auto groupsUpTo = [&](int64_t x, const vector<int64_t> &sum) {   /* groups at boundaries <= x */
        return sum[upper_bound(groupAt.begin(), groupAt.end(), x) - groupAt.begin() - 1];
    };
    auto countLE = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(upper_bound(v.begin(), v.end(), x) - v.begin()); };
    auto countLT = [](const vector<int64_t> &v, int64_t x) { return (int64_t)(lower_bound(v.begin(), v.end(), x) - v.begin()); };
    auto rowsIn = [&](int64_t s, int64_t e) { return countLT(starts, e) - countLE(ends, s); };   /* chains overlapping [s, e) */

    /* ---- blocks.  A block [s, e) holds the reference columns s..e-1 and the insertion
     * groups at boundaries (s, e]: a group where the block is cut ends the left block. */
    size_t nextChain = 0;
    vector<size_t> active;
    vector<BlockRow> block;
    string dna;
    auto emit = [&](int64_t s, int64_t e) {
        while (nextChain < chains.size() && chains[nextChain].start() < e) active.push_back(nextChain++);
        vector<size_t> rows;
        size_t w = 0;
        for (size_t k : active) {
            if (chains[k].end() <= s) continue;
            active[w++] = k;
            rows.push_back(k);
        }
        active.resize(w);
        /* column layout */
        vector<int64_t> refCol((size_t)(e - s));
        map<int64_t, int64_t> groupCol;   /* boundary -> first column */
        int64_t ncol = 0;
        auto gIt = groups.upper_bound(s);
        for (int64_t x = s; x < e; x++) {
            if (x > s && gIt != groups.end() && gIt->first == x) {
                groupCol[x] = ncol;
                ncol += gIt->second.width;
                stats.insertionColumns += gIt->second.width;
                stats.insertedBases += gIt->second.bases;
                ++gIt;
            }
            refCol[(size_t)(x - s)] = ncol++;
        }
        if (gIt != groups.end() && gIt->first == e) {
            groupCol[e] = ncol;
            ncol += gIt->second.width;
            stats.insertionColumns += gIt->second.width;
            stats.insertedBases += gIt->second.bases;
        }
        struct Row {
            const Chain *c;
            int64_t o1, o2;   /* oriented bases [o1, o2) */
        };
        vector<Row> out;
        for (size_t k : rows) {
            const Chain &c = chains[k];
            int64_t o1 = -1, o2 = -1;
            auto span = [&](int64_t a, int64_t b) {
                if (o1 < 0 || a < o1) o1 = a;
                if (b > o2) o2 = b;
            };
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                if (a < b) span(c.ostart[i] + (a - p.srcStart), c.ostart[i] + (b - p.srcStart));
                if (i < c.at.size() && c.at[i] > s && c.at[i] <= e) span(c.oend(i), c.ostart[i + 1]);
            }
            if (o1 >= 0) out.push_back(Row{&c, o1, o2});
        }
        stable_sort(out.begin(), out.end(), [&](const Row &x, const Row &y) { return rowBefore(*x.c, x.o1, *y.c, y.o1); });
        size_t nrows = out.size() + 1;
        if (block.size() < nrows) block.resize(nrows);
        BlockRow &r0 = block[0];
        r0.id = nullptr;
        r0.genome = z.genome;
        r0.seq = z.seq;
        r0.reversed = false;
        r0.start = z.at(s);
        r0.size = z.at(e) - r0.start;
        fetch(z.genome, z.seq, r0.start, r0.start + r0.size, false, dna);
        r0.text.assign((size_t)ncol, '-');
        if (z.all) {
            for (int64_t x = s; x < e; x++) r0.text[(size_t)refCol[(size_t)(x - s)]] = dna[(size_t)(x - s)];
        } else {
            /* the runs from the one holding s, if any */
            size_t i = (size_t)(upper_bound(z.runU.begin(), z.runU.end(), s) - z.runU.begin()), k = 0;
            if (i > 0 && z.runU[i - 1] + z.runLen[i - 1] > s) i--;
            for (; i < z.runU.size() && z.runU[i] < e; i++)
                for (int64_t x = max(s, z.runU[i]); x < min(e, z.runU[i] + z.runLen[i]); x++)
                    r0.text[(size_t)refCol[(size_t)(x - s)]] = dna[k++];
        }
        for (size_t k = 0; k < out.size(); k++) {
            const Row &r = out[k];
            const Chain &c = *r.c;
            BlockRow &br = block[k + 1];
            br.id = &c;
            br.genome = c.genome;
            br.seq = c.seq;
            br.reversed = c.reversed;
            br.start = r.o1;
            br.size = r.o2 - r.o1;
            fetch(c.genome, c.seq, r.o1, r.o2, c.reversed, dna);
            string &text = br.text;
            text.assign((size_t)ncol, '-');
            for (size_t i = 0; i < c.pieces.size(); i++) {
                const LiftPiece &p = c.pieces[i];
                int64_t a = max(s, p.srcStart), b = min(e, p.srcStart + p.length);
                for (int64_t x = a; x < b; x++)
                    text[(size_t)refCol[(size_t)(x - s)]] = dna[(size_t)(c.ostart[i] + (x - p.srcStart) - r.o1)];
                if (i < c.at.size() && c.at[i] > s && c.at[i] <= e) {
                    int64_t base = groupCol[c.at[i]], o = c.oend(i);
                    const vector<int> &cols = c.cols[i];
                    for (size_t j = 0; j < cols.size(); j++) text[(size_t)(base + cols[j])] = dna[(size_t)(o + (int64_t)j - r.o1)];
                }
            }
            stats.alignedBases += r.o2 - r.o1;
        }
        writer.block(block, nrows, ncol);
        stats.blocks++;
        stats.rows += nrows;
        stats.columns += ncol;
    };
    /* close [s, e), splitting at the caps (counted in reference columns) */
    auto close = [&](int64_t s, int64_t e) {
        while (s < e) {
            int64_t w = min(e - s, _opts.maxBlockCols);
            while (w > 1 && (rowsIn(s, s + w) + 1) * w > _opts.maxBlockCells) w = max<int64_t>(1, w / 2);
            emit(s, s + w);
            s += w;
        }
    };

    vector<int64_t> bps;
    for (int64_t x : starts)
        if (x > cs && x < ce) bps.push_back(x);
    for (int64_t x : ends)
        if (x > cs && x < ce) bps.push_back(x);
    for (int64_t x : extra)
        if (x > cs && x < ce) bps.push_back(x);
    sort(bps.begin(), bps.end());
    bps.erase(unique(bps.begin(), bps.end()), bps.end());
    int64_t s = cs;
    for (size_t i = 0; i < bps.size(); i++) {
        int64_t b = bps[i], nb = i + 1 < bps.size() ? bps[i + 1] : ce;
        int64_t endsAt = countLE(ends, b) - countLT(ends, b), startsAt = countLE(starts, b) - countLT(starts, b);
        int64_t activeRows = (countLT(starts, b) - countLT(ends, b)) + startsAt;
        bool cut = activeRows > 0 && (double)(endsAt + startsAt) >= _opts.breakFrac * (double)activeRows;
        if (!cut) {
            /* the block continued to the next breakpoint: its cells and aligned bases,
             * insertion groups included */
            int64_t width = (nb - s) + groupsUpTo(nb, groupCols) - groupsUpTo(s, groupCols);
            int64_t nongap = aligned(nb) - aligned(s) + groupsUpTo(nb, groupBases) - groupsUpTo(s, groupBases);
            int64_t gaps = rowsIn(s, nb) * width - nongap;
            cut = (double)gaps > _opts.gapBudget * (double)nongap;
        }
        if (cut) {
            close(s, b);
            s = b;
        }
    }
    close(s, ce);
}

void addStats(MafStats &a, const MafStats &b) {
    a.blocks += b.blocks;
    a.rows += b.rows;
    a.columns += b.columns;
    a.alignedBases += b.alignedBases;
    a.insertionColumns += b.insertionColumns;
    a.insertedBases += b.insertedBases;
    a.joinsRejected += b.joinsRejected;
    a.backboneColumns += b.backboneColumns;
    a.placedColumns += b.placedColumns;
}

bool endsWith(const string &s, const string &t) { return s.size() >= t.size() && !s.compare(s.size() - t.size(), t.size(), t); }

/* --universalBed: the axis in order, one BED6 line per maximal run of one sequence --
 * genome.sequence, start, end, first column, depth, strand -- a run cut by a chunk edge
 * written once.  BGZF when the name ends in .gz. */
class AxisBed {
  public:
    AxisBed(const string &path, const UniversalIndex &index)
        : _file(path), _bgzf(endsWith(path, ".gz")), _out(_file, _bgzf), _idx(index) {}
    void add(const vector<Segment> &segs, int64_t u0) {
        for (const Segment &sg : segs) {
            int64_t u = u0 + sg.u;
            if (_have && sg.g == _run.g && sg.seq == _run.seq && sg.rev == _run.rev && sg.depth == _run.depth && u == _next &&
                (sg.rev ? sg.pos + sg.len == _run.pos : sg.pos == _run.pos + _run.len)) {
                if (sg.rev) _run.pos = sg.pos;
                _run.len += sg.len;
                _next += sg.len;
                continue;
            }
            flush();
            _have = true;
            _run = sg;
            _run.u = u;
            _next = u + sg.len;
        }
    }
    void finish() {
        flush();
        _out.put(_line);
        _out.finish();
        if (_bgzf) writeBgzfEof(_file);
        _file.close();
    }

  private:
    FileSink _file;
    bool _bgzf;
    TextOut _out;
    const UniversalIndex &_idx;
    bool _have = false;
    Segment _run;
    int64_t _next = 0;
    string _line;
    void flush() {
        if (!_have) return;
        const Genome *g = _idx.genome(_run.g);
        _line += rowName(g, _run.seq, false) + "\t" + to_string(_run.pos) + "\t" + to_string(_run.pos + _run.len) + "\t" +
                 to_string(_run.u) + "\t" + to_string(_run.depth) + (_run.rev ? "\t-\n" : "\t+\n");
        _have = false;
        if (_line.size() >= ((size_t)1 << 16)) {
            _out.put(_line);
            _line.clear();
        }
    }
};

/* A whole universal export wrote every base due exactly once: per genome, the moments
 * of what its rows held against what they should hold -- the backbone all of it; an
 * ancestor its top intervals (all with allAncestorRows); a leaf all but its top
 * intervals (all with leafInsertions); an excluded leaf nothing. */
void checkOnce(const UniversalIndex &idx, const MafOptions &opts, const map<const Genome *, Moments> &got) {
    set<string> excluded(opts.exclude.begin(), opts.exclude.end());
    for (int r = 0; r < idx.genomes(); r++) {
        const Genome *g = idx.genome(r);
        const bool leaf = g->children().empty();
        Moments want;
        if (!excluded.count(g->name())) {
            if (r == 0 || (leaf ? opts.leafInsertions : opts.allAncestorRows)) {
                want = idx.allBases(r);
            } else if (!leaf) {
                want = idx.topBases(r);
            } else {
                want = idx.allBases(r);
                want.subtract(idx.topBases(r));
            }
        }
        auto h = got.find(g);
        const Moments have = h == got.end() ? Moments() : h->second;
        if (!(have == want))
            throw runtime_error("internal error: --universal wrote " + to_string(have.n) + " bases of " + g->name() + " where " +
                                to_string(want.n) + " were due, each once");
    }
}
}

/* ---------------------------------------------------------------- the driver */

MafStats exportMaf(const Alignment &aln, const Genome *ref, const MafOptions &opts, const MafFormat &format,
                   const vector<RefInterval> &intervals, const string &path, int threads) {
    { TreeWalk check(aln, ref, walkOptions(opts)); }   /* bad options fail before any output */
    const map<const Genome *, int> ranks = rowRanks(aln, ref);
    /* chunks: of the reference intervals; universal, of the axis the planner lays out */
    vector<RefInterval> chunks;
    unique_ptr<UniversalIndex> index;
    unique_ptr<UniversalPlanner> plan;
    if (opts.universal) {
        set<const Genome *> excluded;
        for (const string &x : opts.exclude) excluded.insert(aln.genome(x));
        index.reset(new UniversalIndex(aln, ref, opts.leafInsertions, opts.dupes, excluded, opts.universalWhole, threads));
        vector<TopRange> ranges;
        for (const RefInterval &iv : intervals) ranges.push_back(TopRange{0, iv.seq, iv.start, iv.end});
        plan.reset(new UniversalPlanner(*index, ranges, opts.universalWhole, opts.chunk));
        plan->run();
    } else {
        chunks = chunkIntervals(intervals, opts.chunk);
    }
    const size_t nChunks = plan ? plan->chunks() : chunks.size();
    unique_ptr<AxisBed> bed;
    if (plan && !opts.universalBed.empty()) bed.reset(new AxisBed(opts.universalBed, *index));
    FileSink file(path);
    {
        TextOut header(file, format.bgzf);
        /* universal: what the axis is, as tags; and the tree, as hal2maf gives it (taffy
         * skips # lines; its universal index needs it) */
        string tags;
        if (opts.universal) {
            const char *eq = format.taf ? ":" : "=";
            tags = string(" universal") + eq + "flank1 leafInsertions" + eq + (opts.leafInsertions ? "1" : "0") + " ancestorRows" + eq +
                   (opts.allAncestorRows ? "all" : "own") + " noDupes" + eq + (opts.dupes ? "0" : "1") + " chunk" + eq +
                   to_string(opts.chunk);
            if (ref != aln.root()) tags += string(" rootGenome") + eq + ref->name();
        }
        if (format.taf)
            header.put(string(format.runLength ? "#taf run_length_encode_bases:1 version:1 scoring:N/A" : "#taf version:1 scoring:N/A") +
                       tags + "\n# hal " + aln.manifest().newick + "\n");
        else header.put("##maf version=1 scoring=N/A" + tags + "\n# hal " + aln.manifest().newick + "\n\n");
        header.finish();
    }
    MafStats total;
    map<const Genome *, Moments> moments;
    size_t prevRows = 0;   /* TAF: rows the last chunk written ended with */
    auto firstLine = [&](const ChunkResult &r) {   /* its own BGZF block(s), in either path */
        TextOut first(file, format.bgzf);
        first.put(tafFirstLine(r, prevRows));
        first.finish();
    };
    /* one exporter per thread, made when the thread first needs it */
    vector<unique_ptr<ChunkExporter>> exporters((size_t)max(1, threads));
    /* chunk k by thread w to out; direct: out is the output, and the chunk's first line
     * goes there as soon as it is known */
    auto exportChunk = [&](size_t k, int w, ByteSink &out, ChunkResult &r, bool direct) {
        TextOut body(out, format.bgzf);
        unique_ptr<ChunkExporter> &e = exporters[(size_t)w];
        if (!e) e.reset(new ChunkExporter(aln, ref, opts, ranks));
        ChunkExporter &ex = *e;
        auto run = [&](BlockWriter &bw) {
            if (!plan) {
                ex.run(chunks[k], bw);
                return;
            }
            plan->segments(k, r.axis);
            const UniversalPlanner::Chunk &c = plan->chunk(k);
            const TopRange &t = plan->ranges()[c.range];
            ex.runUniversal(UniChunk{index.get(), &r.axis, index->genome(t.g), t.seq, c.rangePos, c.cols}, bw);
            r.moments.swap(ex.moments);
        };
        if (format.taf) {
            TafBlockWriter w(body, format.runLength, opts.onlySequenceNames, opts.universal, r);
            if (direct) w.onFirstLine([&]() { firstLine(r); });
            run(w);
            w.finish();
        } else {
            MafBlockWriter w(body, opts.onlySequenceNames);
            run(w);
        }
        body.finish();
        r.stats = ex.stats;
    };
    /* what is kept of a chunk once it is written */
    auto account = [&](size_t k, ChunkResult &r) {
        prevRows = r.lastRows;
        addStats(total, r.stats);
        if (bed) bed->add(r.axis, plan->chunk(k).u0);
        for (auto &m : r.moments) moments[m.first].add(m.second);
    };
    if (threads <= 1 || nChunks <= 1) {
        for (size_t k = 0; k < nChunks; k++) {
            ChunkResult r;
            exportChunk(k, 0, file, r, true);
            account(k, r);
        }
    } else {
        /* chunks wait for the ones before them beside the output, or where temporary
         * files go when none can be made there (stdout, /dev/stdout, >(...)) */
        vector<string> tmpDirs;
        size_t k = path.rfind('/');
        if (path != "stdout") tmpDirs.push_back(k == string::npos ? "." : k == 0 ? "/" : path.substr(0, k));
        const char *t = getenv("TMPDIR");
        tmpDirs.push_back(t && *t ? t : "/tmp");
        orderedParallel<ChunkResult>(
            nChunks, threads,
            [&](size_t k, int w) {
                ChunkResult r;
                r.body.reset(new ChunkBuffer(tmpDirs));
                exportChunk(k, w, *r.body, r, false);
                return r;
            },
            [&](size_t k, ChunkResult &r) {
                if (format.taf) firstLine(r);
                r.body->copyTo(file);
                r.body.reset();
                account(k, r);
            });
    }
    if (format.bgzf) writeBgzfEof(file);
    file.close();
    if (bed) bed->finish();
    if (plan) {
        total.placedIntervals = plan->visited();
        if (opts.universalWhole) {
            total.extraSequences = index->flankless().size();
            checkOnce(*index, opts, moments);
            total.onceChecked = true;
        }
    }
    return total;
}
}
