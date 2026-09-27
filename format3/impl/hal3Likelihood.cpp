/*
 * HAL format 3: per-base likelihoods, the .hlk file.  See hal3Likelihood.h.
 */
#include "hal3Likelihood.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <unistd.h>
#include <unordered_map>
#include <zlib.h>

using namespace std;

namespace hal3 {

namespace {

const char MAGIC[8] = {'\x89', 'H', 'L', 'K', '\r', '\n', '\x1a', '\n'};
const uint32_t VERSION = 1;
const size_t HEADER = 64, DIRENT = 40;

void put32(string &s, uint32_t v) {
    for (int k = 0; k < 4; k++) s += (char)(v >> (8 * k));
}
void put64(string &s, uint64_t v) {
    for (int k = 0; k < 8; k++) s += (char)(v >> (8 * k));
}
uint32_t get32(const unsigned char *p) { return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24; }
uint64_t get64(const unsigned char *p) { return (uint64_t)get32(p) | (uint64_t)get32(p + 4) << 32; }
uint32_t crc(const void *p, size_t n) { return (uint32_t)crc32(crc32(0L, Z_NULL, 0), (const Bytef *)p, (uInt)n); }

uint32_t key(const LikCode &c) { return (uint32_t)c.d | (uint32_t)c.k0 << 8 | (uint32_t)c.k1 << 16 | (uint32_t)c.k2 << 24; }

/* Bounds-checked reading of a byte buffer. */
struct Cursor {
    const unsigned char *p, *end;
    const string &what;
    void need(size_t n) const {
        if ((size_t)(end - p) < n) throw runtime_error(what + ": truncated footer");
    }
    uint32_t u32() {
        need(4);
        uint32_t v = get32(p);
        p += 4;
        return v;
    }
    uint64_t u64() {
        need(8);
        uint64_t v = get64(p);
        p += 8;
        return v;
    }
    string str(size_t maxLen) {
        uint32_t n = u32();
        if (n > maxLen) throw runtime_error(what + ": a footer string is too long");
        need(n);
        string s((const char *)p, n);
        p += n;
        return s;
    }
};

string deflateRaw(const string &in, uint8_t &coding) {
    coding = 0;
    if (in.empty()) return in;
    z_stream z;
    memset(&z, 0, sizeof z);
    if (deflateInit2(&z, 9, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) throw runtime_error("deflateInit2 failed");
    string out(deflateBound(&z, (uLong)in.size()), '\0');
    z.next_in = (Bytef *)in.data();
    z.avail_in = (uInt)in.size();
    z.next_out = (Bytef *)&out[0];
    z.avail_out = (uInt)out.size();
    int rc = deflate(&z, Z_FINISH);
    size_t n = z.total_out;
    deflateEnd(&z);
    if (rc != Z_STREAM_END) throw runtime_error("deflate failed");
    if (n >= in.size()) return in;
    out.resize(n);
    coding = 1;
    return out;
}

void inflateRaw(const unsigned char *p, size_t n, size_t rawLen, string &out, const string &what) {
    out.assign(rawLen, '\0');
    if (rawLen == 0) {
        if (n != 0) throw runtime_error(what + ": a stream inflates to nothing");
        return;
    }
    z_stream z;
    memset(&z, 0, sizeof z);
    if (inflateInit2(&z, -15) != Z_OK) throw runtime_error("inflateInit2 failed");
    z.next_in = (Bytef *)p;
    z.avail_in = (uInt)n;
    z.next_out = (Bytef *)&out[0];
    z.avail_out = (uInt)rawLen;
    int rc = inflate(&z, Z_FINISH);
    size_t got = z.total_out;
    inflateEnd(&z);
    if (rc != Z_STREAM_END || got != rawLen) throw runtime_error(what + ": a stream does not inflate to its size");
}

void uvarint(string &s, uint64_t v) {
    while (v >= 0x80) {
        s += (char)(v | 0x80);
        v >>= 7;
    }
    s += (char)v;
}

/* 10^(-k/q) for k in 0..254, and 0 for 255 */
struct Powers {
    double v[256];
    explicit Powers(int q) {
        for (int k = 0; k < 255; k++) v[k] = pow(10.0, -(double)k / q);
        v[255] = 0;
    }
};

void decodeWith(const LikCode &c, int r, const Powers &pw, double p[4]) {
    int m = (r < 4 ? r : 0) ^ (c.d & 3);
    p[m] = 1;
    p[m ^ 2] = pw.v[c.k0];
    p[m ^ 1] = pw.v[c.k1];
    p[m ^ 3] = pw.v[c.k2];
    double s = p[0] + p[1] + p[2] + p[3];
    for (int b = 0; b < 4; b++) p[b] /= s;
}
}

int likBase(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

LikCode encodeLikelihood(const double p[4], int r, int q) {
    int rr = r < 4 ? r : 0;
    double mx = -1;
    for (int b = 0; b < 4; b++) {
        if (!(p[b] >= 0) || !isfinite(p[b])) return LikCode{0, 0, 0, 0};   /* nonsense: no information */
        mx = max(mx, p[b]);
    }
    if (!(mx > 0)) return LikCode{0, 0, 0, 0};
    int m = rr;
    if (p[rr] != mx)
        for (m = 0; p[m] != mx; m++) {
        }
    uint8_t k[3];
    const int xs[3] = {m ^ 2, m ^ 1, m ^ 3};
    for (int j = 0; j < 3; j++) {
        double px = p[xs[j]];
        if (!(px > 0)) {
            k[j] = 255;
            continue;
        }
        double v = q * log10(mx / px);
        k[j] = v >= 253.5 ? 254 : (uint8_t)lround(v);   /* nonzero: never stored as 0 */
    }
    return LikCode{(uint8_t)(m ^ rr), k[0], k[1], k[2]};
}

void decodeLikelihood(const LikCode &c, int r, int q, double p[4]) { decodeWith(c, r, Powers(q), p); }

/* ---------------------------------------------------------------- reading */

LikFile::~LikFile() {
    if (_ownIo && _io) lmap_io_close(_io);
}

unique_ptr<LikFile> LikFile::open(lmap_io *io, bool ownIo, const string &what) {
    unique_ptr<LikFile> f(new LikFile());
    f->_io = io;
    f->_ownIo = ownIo;
    f->_what = what;
    int64_t size = lmap_io_size(io);
    unsigned char h[HEADER];
    if (size < (int64_t)HEADER || lmap_pread(io, h, 0, HEADER) != 0) throw runtime_error(what + ": cannot read the header");
    if (memcmp(h, MAGIC, 8) != 0) throw runtime_error(what + " is not a likelihood (.hlk) file");
    if (get32(h + 60) != crc(h, 60)) throw runtime_error(what + ": the header is damaged");
    if (get32(h + 8) != VERSION) throw runtime_error(what + ": version " + to_string(get32(h + 8)) + " is not supported");
    if (get32(h + 12) != 0) throw runtime_error(what + ": uses features this reader does not know");
    uint64_t fo = get64(h + 16), fl = get64(h + 24);
    if (fo < HEADER || fl > (uint64_t)size || fo > (uint64_t)size - fl) throw runtime_error(what + ": the footer is outside the file");
    vector<unsigned char> foot((size_t)fl);
    if (fl && lmap_pread(io, foot.data(), (int64_t)fo, (int64_t)fl) != 0) throw runtime_error(what + ": cannot read the footer");
    if (crc(foot.data(), foot.size()) != get32(h + 32)) throw runtime_error(what + ": the footer is damaged");
    Cursor c{foot.data(), foot.data() + foot.size(), f->_what};
    f->_q = (int)c.u32();
    f->_chunkBases = c.u32();
    uint32_t nd = c.u32();
    if (f->_q < 1 || f->_q > LIK_MAX_STEPS || f->_chunkBases == 0 || f->_chunkBases > LIK_MAX_CHUNK || nd < 1 || nd > 255)
        throw runtime_error(what + ": bad parameters");
    for (uint32_t i = 0; i < nd; i++) {
        c.need(4);
        LikCode k{c.p[0], c.p[1], c.p[2], c.p[3]};
        c.p += 4;
        if (k.d > 3) throw runtime_error(what + ": a bad dictionary code");
        f->_dict.push_back(k);
    }
    for (uint32_t i = 0, n = c.u32(); i < n; i++) {
        string k = c.str(1 << 20);
        f->_meta[k] = c.str(1 << 30);
    }
    uint64_t chunks = 0;
    for (uint32_t i = 0, n = c.u32(); i < n; i++) {
        Seq s;
        s.name = c.str(1 << 16);
        s.length = (int64_t)c.u64();
        s.firstChunk = c.u64();
        if (s.length < 0 || s.firstChunk != chunks) throw runtime_error(what + ": a bad sequence table");
        chunks += ((uint64_t)s.length + f->_chunkBases - 1) / f->_chunkBases;
        f->_seqs.push_back(s);
    }
    f->_nChunks = c.u64();
    f->_dirOffset = c.u64();
    f->_dirCrc = c.u32();
    if (f->_nChunks != chunks) throw runtime_error(what + ": the chunk count does not match the sequences");
    if (f->_dirOffset < HEADER || f->_dirOffset > (uint64_t)size || f->_nChunks > ((uint64_t)size - f->_dirOffset) / DIRENT)
        throw runtime_error(what + ": the directory is outside the file");
    if (c.p != c.end) throw runtime_error(what + ": bytes after the footer's end");
    return f;
}

string LikFile::meta(const string &key) const {
    auto i = _meta.find(key);
    return i == _meta.end() ? "" : i->second;
}

void LikFile::decodeChunk(uint64_t chunk, int64_t bases, vector<LikCode> &out) const {
    string where = _what + ", chunk " + to_string(chunk);
    unsigned char e[DIRENT];
    if (lmap_pread(_io, e, (int64_t)(_dirOffset + chunk * DIRENT), DIRENT) != 0) throw runtime_error(where + ": cannot read its entry");
    uint64_t off = get64(e);
    uint32_t listed = get32(e + 8), gz = get32(e + 12), iz = get32(e + 16), ez = get32(e + 20);
    uint32_t graw = get32(e + 24), eraw = get32(e + 28);
    uint8_t cg = e[32], ci = e[33], ce = e[34];
    uint32_t want = get32(e + 36);
    uint64_t n = (uint64_t)gz + iz + ez;
    /* every size bounded by the chunk's bases before anything is allocated: a gap is a
     * uvarint of at most 4 bytes (chunks hold at most 2^24 bases), a code 4 bytes, and a
     * deflated stream is smaller than its raw bytes (the writer stores it otherwise) */
    auto sized = [](uint8_t coding, uint32_t stored, uint64_t raw) { return coding ? stored < raw : stored == raw; };
    if (cg > 1 || ci > 1 || ce > 1 || (int64_t)listed > bases || graw > 4 * (uint64_t)listed || graw < listed ||
        eraw > 4 * (uint64_t)listed || eraw % 4 || !sized(cg, gz, graw) || !sized(ci, iz, listed) || !sized(ce, ez, eraw) ||
        off < HEADER || n > _dirOffset || off > _dirOffset - n)
        throw runtime_error(where + ": a bad entry");
    vector<unsigned char> raw((size_t)n);
    if (n && lmap_pread(_io, raw.data(), (int64_t)off, (int64_t)n) != 0) throw runtime_error(where + ": cannot read it");
    if (crc(raw.data(), raw.size()) != want) throw runtime_error(where + ": damaged (CRC)");
    string G, I, E;
    const unsigned char *p = raw.data();
    auto stream = [&](uint8_t coding, uint32_t stored, size_t rawLen, string &s) {
        if (coding) inflateRaw(p, stored, rawLen, s, where);
        else if (stored != rawLen) throw runtime_error(where + ": a stored stream has the wrong size");
        else s.assign((const char *)p, stored);
        p += stored;
    };
    stream(cg, gz, graw, G);
    stream(ci, iz, listed, I);
    stream(ce, ez, eraw, E);
    out.assign((size_t)bases, _dict[0]);
    size_t gi = 0, ei = 0;
    int64_t pos = -1;
    for (uint32_t i = 0; i < listed; i++) {
        uint64_t g = 0;
        for (int shift = 0;; shift += 7) {
            if (gi >= G.size() || shift > 63) throw runtime_error(where + ": a bad gap stream");
            unsigned char b = (unsigned char)G[gi++];
            g |= (uint64_t)(b & 0x7f) << shift;
            if (!(b & 0x80)) break;
        }
        if (g >= (uint64_t)bases) throw runtime_error(where + ": a gap past the chunk");
        pos += 1 + (int64_t)g;
        if (pos >= bases) throw runtime_error(where + ": a base past the chunk");
        unsigned char idx = (unsigned char)I[i];
        if (idx == 255) {
            if (ei + 4 > E.size()) throw runtime_error(where + ": too few codes in E");
            LikCode k{(uint8_t)E[ei], (uint8_t)E[ei + 1], (uint8_t)E[ei + 2], (uint8_t)E[ei + 3]};
            if (k.d > 3) throw runtime_error(where + ": a bad code");
            out[(size_t)pos] = k;
            ei += 4;
        } else {
            if (idx >= _dict.size()) throw runtime_error(where + ": a dictionary index past the dictionary");
            out[(size_t)pos] = _dict[idx];
        }
    }
    if (gi != G.size() || ei != E.size()) throw runtime_error(where + ": streams longer than their bases");
}

void LikFile::codes(uint32_t seq, int64_t start, int64_t end, vector<LikCode> &out, LikCache *cache) const {
    if (seq >= _seqs.size() || start < 0 || end < start || end > _seqs[seq].length)
        throw runtime_error(_what + ": no bases " + to_string(start) + "-" + to_string(end) + " in sequence " + to_string(seq));
    out.resize((size_t)(end - start));
    const Seq &s = _seqs[seq];
    vector<LikCode> local;
    for (int64_t x = start; x < end;) {
        int64_t c = x / _chunkBases, cs = c * _chunkBases, ce = min<int64_t>(s.length, cs + _chunkBases);
        uint64_t global = s.firstChunk + (uint64_t)c;
        const vector<LikCode> *src;
        if (cache && cache->file == this && cache->chunk == (int64_t)global) {
            src = &cache->codes;
        } else if (cache) {
            cache->file = nullptr;
            decodeChunk(global, ce - cs, cache->codes);
            cache->file = this;
            cache->chunk = (int64_t)global;
            src = &cache->codes;
        } else {
            decodeChunk(global, ce - cs, local);
            src = &local;
        }
        int64_t y = min(end, ce);
        copy(src->begin() + (x - cs), src->begin() + (y - cs), out.begin() + (x - start));
        x = y;
    }
}

void LikFile::probabilities(uint32_t seq, int64_t start, int64_t end, const char *dna, bool reversed, double *out,
                            LikCache *cache) const {
    vector<LikCode> cs;
    codes(seq, start, end, cs, cache);
    Powers pw(_q);
    size_t n = cs.size();
    for (size_t i = 0; i < n; i++) {
        size_t j = reversed ? n - 1 - i : i;   /* the forward base */
        int r = likBase(dna[j]);
        if (reversed) r = (r < 4 ? r : 0) ^ 3;   /* an N counts as A, so its complement T */
        decodeWith(cs[j], r, pw, out + 4 * i);
    }
}

vector<string> LikFile::verify() const {
    vector<string> problems;
    uint64_t dirBytes = _nChunks * DIRENT;
    vector<unsigned char> dir((size_t)dirBytes);
    if (dirBytes && lmap_pread(_io, dir.data(), (int64_t)_dirOffset, (int64_t)dirBytes) != 0) {
        problems.push_back(_what + ": cannot read the directory");
        return problems;
    }
    if (crc(dir.data(), dir.size()) != _dirCrc) problems.push_back(_what + ": the directory is damaged (CRC)");
    for (size_t i = 0; i < _dict.size(); i++)
        for (size_t j = 0; j < i; j++)
            if (_dict[i] == _dict[j]) problems.push_back(_what + ": dictionary entries " + to_string(j) + " and " + to_string(i) + " are the same");
    vector<LikCode> codes;
    for (const Seq &s : _seqs)
        for (int64_t c = 0; c * (int64_t)_chunkBases < s.length; c++) {
            try {
                decodeChunk(s.firstChunk + (uint64_t)c, min<int64_t>(_chunkBases, s.length - c * (int64_t)_chunkBases), codes);
            } catch (exception &e) {
                problems.push_back(e.what());
                if (problems.size() > 20) return problems;
            }
        }
    return problems;
}

/* ---------------------------------------------------------------- writing */

LikDictionary LikDictionary::choose(const LikCode *codes, size_t n) {
    unordered_map<uint32_t, uint64_t> count;
    for (size_t i = 0; i < n; i++) count[key(codes[i])]++;
    vector<pair<uint64_t, uint32_t>> v;
    for (auto &kc : count) v.push_back({kc.second, kc.first});
    sort(v.begin(), v.end(), [](const pair<uint64_t, uint32_t> &a, const pair<uint64_t, uint32_t> &b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;
    });
    if (v.size() > 255) v.resize(255);
    if (v.empty()) v.push_back({0, key(LikCode{0, 255, 255, 255})});
    LikDictionary d;
    for (size_t i = 0; i < v.size(); i++) {
        uint32_t k = v[i].second;
        d.codes.push_back(LikCode{(uint8_t)k, (uint8_t)(k >> 8), (uint8_t)(k >> 16), (uint8_t)(k >> 24)});
        d.index[k] = (uint8_t)i;
    }
    return d;
}

LikChunk encodeLikChunk(const LikCode *codes, size_t n, const LikDictionary &dict) {
    string G, I, E;
    LikChunk c;
    c.bases = (uint32_t)n;
    int64_t last = -1;
    for (size_t i = 0; i < n; i++) {
        const LikCode &x = codes[i];
        if (x == dict.codes[0]) continue;
        uvarint(G, (uint64_t)((int64_t)i - last - 1));
        last = (int64_t)i;
        c.listed++;
        auto k = dict.index.find(key(x));
        if (k != dict.index.end()) {
            I += (char)k->second;
        } else {
            I += (char)255;
            E += (char)x.d;
            E += (char)x.k0;
            E += (char)x.k1;
            E += (char)x.k2;
        }
    }
    string g = deflateRaw(G, c.cg), in = deflateRaw(I, c.ci), e = deflateRaw(E, c.ce);
    c.gz = (uint32_t)g.size();
    c.iz = (uint32_t)in.size();
    c.ez = (uint32_t)e.size();
    c.graw = (uint32_t)G.size();
    c.eraw = (uint32_t)E.size();
    c.bytes = g + in + e;
    c.crc = crc(c.bytes.data(), c.bytes.size());
    return c;
}

LikWriter::LikWriter(const string &path, const vector<pair<string, int64_t>> &seqs, const map<string, string> &meta,
                     const LikDictionary &dict, int q, uint32_t chunkBases)
    : _path(path), _partial(path + ".partial"), _seqs(seqs), _meta(meta), _dict(dict), _q(q), _chunkBases(chunkBases) {
    if (q < 1 || q > LIK_MAX_STEPS || chunkBases == 0 || chunkBases > LIK_MAX_CHUNK || dict.codes.empty() || dict.codes.size() > 255)
        throw runtime_error(path + ": bad likelihood file parameters");
    uint64_t chunks = 0;
    for (auto &s : _seqs) {
        if (s.second < 0) throw runtime_error(path + ": a negative sequence length");
        _firstChunk.push_back(chunks);
        chunks += ((uint64_t)s.second + chunkBases - 1) / chunkBases;
    }
    _f = fopen(_partial.c_str(), "wb");
    if (!_f) throw runtime_error("cannot write " + _partial + ": " + strerror(errno));
    char zero[HEADER] = {0};
    put(zero, HEADER);
    while (_seq < _seqs.size() && _seqs[_seq].second == 0) _seq++;
}

LikWriter::~LikWriter() {
    if (_f) {
        fclose(_f);
        unlink(_partial.c_str());
    }
}

void LikWriter::put(const void *p, size_t n) {
    if (n && fwrite(p, 1, n, _f) != n) throw runtime_error("cannot write " + _partial);
}

int64_t LikWriter::nextChunkBases() const {
    if (_seq >= _seqs.size()) return 0;
    return min<int64_t>(_chunkBases, _seqs[_seq].second - _inSeq);
}

void LikWriter::addChunk(const LikChunk &c) {
    if (!_f) throw runtime_error(_path + ": written to after close");
    if (!_chunk.empty() || (int64_t)c.bases != nextChunkBases() || c.bases == 0)
        throw runtime_error(_path + ": a chunk of " + to_string(c.bases) + " bases where " + to_string(nextChunkBases()) + " belong");
    put(c.bytes.data(), c.bytes.size());
    put64(_dir, _offset);
    put32(_dir, c.listed);
    put32(_dir, c.gz);
    put32(_dir, c.iz);
    put32(_dir, c.ez);
    put32(_dir, c.graw);
    put32(_dir, c.eraw);
    _dir += (char)c.cg;
    _dir += (char)c.ci;
    _dir += (char)c.ce;
    _dir += (char)0;
    put32(_dir, c.crc);
    _offset += c.bytes.size();
    _nChunks++;
    _inSeq += c.bases;
    if (_inSeq == _seqs[_seq].second) {
        _inSeq = 0;
        for (_seq++; _seq < _seqs.size() && _seqs[_seq].second == 0; _seq++) {
        }
    }
}

void LikWriter::add(const LikCode *codes, size_t n) {
    for (size_t i = 0; i < n; i++) {
        if (nextChunkBases() == 0) throw runtime_error(_path + ": more bases than the sequences hold");
        _chunk.push_back(codes[i]);
        if ((int64_t)_chunk.size() == nextChunkBases()) {
            vector<LikCode> full;
            full.swap(_chunk);
            addChunk(encodeLikChunk(full.data(), full.size(), _dict));
        }
    }
}

void LikWriter::close() {
    if (!_f) throw runtime_error(_path + ": closed twice");
    if (_seq < _seqs.size() || !_chunk.empty()) throw runtime_error(_path + ": not every base was added");
    uint64_t dirOffset = _offset;
    put(_dir.data(), _dir.size());
    _offset += _dir.size();
    string foot;
    put32(foot, (uint32_t)_q);
    put32(foot, _chunkBases);
    put32(foot, (uint32_t)_dict.codes.size());
    for (const LikCode &c : _dict.codes) {
        foot += (char)c.d;
        foot += (char)c.k0;
        foot += (char)c.k1;
        foot += (char)c.k2;
    }
    put32(foot, (uint32_t)_meta.size());
    for (auto &kv : _meta) {
        put32(foot, (uint32_t)kv.first.size());
        foot += kv.first;
        put32(foot, (uint32_t)kv.second.size());
        foot += kv.second;
    }
    put32(foot, (uint32_t)_seqs.size());
    for (size_t i = 0; i < _seqs.size(); i++) {
        put32(foot, (uint32_t)_seqs[i].first.size());
        foot += _seqs[i].first;
        put64(foot, (uint64_t)_seqs[i].second);
        put64(foot, _firstChunk[i]);
    }
    put64(foot, _nChunks);
    put64(foot, dirOffset);
    put32(foot, crc(_dir.data(), _dir.size()));
    uint64_t footOffset = _offset;
    put(foot.data(), foot.size());
    string h(MAGIC, 8);
    put32(h, VERSION);
    put32(h, 0);
    put64(h, footOffset);
    put64(h, foot.size());
    put32(h, crc(foot.data(), foot.size()));
    h.resize(60, '\0');
    put32(h, crc(h.data(), 60));
    if (fseek(_f, 0, SEEK_SET) != 0) throw runtime_error("cannot write " + _partial);
    put(h.data(), h.size());
    FILE *f = _f;
    _f = nullptr;
    if (fclose(f) != 0) {
        unlink(_partial.c_str());
        throw runtime_error("cannot write " + _partial);
    }
    if (rename(_partial.c_str(), _path.c_str()) != 0) {
        unlink(_partial.c_str());
        throw runtime_error("cannot rename " + _partial + " to " + _path + ": " + strerror(errno));
    }
}
}
