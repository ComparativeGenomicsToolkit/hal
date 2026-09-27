/*
 * HAL format 3: per-base likelihoods -- a genome's "lik" member, a .hlk file.
 *
 * For every base of a genome, the probabilities of A, C, G and T.  What they mean is
 * recorded in the file (meta "kind"):
 *   subtree    the likelihood of the data below the genome given each base, normalized
 *              per base: what building an alignment bottom-up propagates to the parent
 *              (Felsenstein's pruning), and, under a uniform prior, the probability of
 *              each base given the genome's subtree;
 *   posterior  the probability of each base given all the data (after a top-down pass);
 *   anything else, as its producer defines it.
 * A genome with no lik member is certain: 1 for its .2bit base, 0 for the others (and
 * for N, no information: 1/4 each).
 *
 * The code of a base.  With bases numbered A0 C1 G2 T3 and r the genome's .2bit base (N
 * or any other letter counts as A), the most probable base is m = r XOR d, and the other
 * three, in the order m^2 (the transition partner), m^1, m^3, each get
 *   k = min(254, round(q log10(p_m / p_x))),   or 255 for probability exactly 0,
 * q being the file's steps per decade (default 10: 0.1 decade, or 1 dB, a step; at most
 * 25, so a file spans at least 10 decades -- a base less likely than that is stored as
 * that).  Where several bases tie for the most probable, m is r if it is one of them,
 * else the first.  Complementing a base is x^3, so a code reads the same on either
 * strand: the reverse strand's vector is decoded from the same code with the complement
 * of r (of A, for an N).
 *
 * The file, little-endian:
 *   header, 64 bytes: magic "\x89HLK\r\n\x1a\n", u32 version 1, u32 flags (0; a reader
 *     refuses bits it does not know), u64 footer offset, u64 footer length, u32 footer
 *     CRC-32, zero padding, and at 60 the CRC-32 of bytes [0, 60).
 *   the chunks, then their directory, then the footer.
 *   footer: u32 q; u32 chunk bases; u32 n (1..255) and the dictionary, n codes of 4
 *     bytes (d, k0, k1, k2); u32 meta count, then per entry u32 key length, key, u32
 *     value length, value; u32 sequence count, then per sequence u32 name length, name,
 *     u64 length, u64 first chunk; u64 chunk count, u64 directory offset, u32 directory
 *     CRC-32.
 *   directory: per chunk 40 bytes: u64 offset, u32 bases listed, u32 stored bytes of G,
 *     of I and of E, u32 bytes of G and of E once inflated, u8 G, I and E coding (0
 *     stored, 1 raw deflate), u8 0, u32 CRC-32 of the chunk's stored bytes.
 * A chunk holds `chunk bases` (at most 2^24) bases of one sequence (a sequence's last
 * chunk fewer), in three streams back to back.  A base whose code is dictionary entry 0 -- the genome's
 * commonest -- is not listed; the others are:
 *   G  uvarint gaps: the first listed base is g after the chunk's start, each next one
 *      1 + g after the one before;
 *   I  a byte per listed base: its dictionary index, or 255 for a code in E;
 *   E  4 bytes per code not in the dictionary.
 * Reads go through lmap_io (local, in a pack, over HTTP), a chunk at a time.
 */
#ifndef _HAL3LIKELIHOOD_H
#define _HAL3LIKELIHOOD_H

#include "lmap_io.h"

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace hal3 {

struct LikCode {
    uint8_t d, k0, k1, k2;
    bool operator==(const LikCode &o) const { return d == o.d && k0 == o.k0 && k1 == o.k1 && k2 == o.k2; }
    bool operator!=(const LikCode &o) const { return !(*this == o); }
};

/* A base's number, A0 C1 G2 T3 in either case, or 4 for anything else. */
int likBase(char c);

/* Encode or decode one base's vector p (A, C, G, T; decode writes it normalized), r
 * being the .2bit base's likBase(). */
LikCode encodeLikelihood(const double p[4], int r, int q);
void decodeLikelihood(const LikCode &c, int r, int q, double p[4]);

const int LIK_MAX_STEPS = 25;
const uint32_t LIK_MAX_CHUNK = (uint32_t)1 << 24;

/* The chunk a reader decoded last, to reuse for nearby reads.  One per thread (it
 * remembers its file). */
struct LikCache {
    const void *file = nullptr;
    int64_t chunk = -1;
    std::vector<LikCode> codes;
};

/* A file's dictionary: up to 255 codes, the commonest first. */
struct LikDictionary {
    std::vector<LikCode> codes;
    std::map<uint32_t, uint8_t> index;
    /* the commonest codes of these (the commonest first; ties by code) */
    static LikDictionary choose(const LikCode *codes, size_t n);
};

/* One chunk, encoded: its stored bytes and directory fields.  encodeLikChunk depends on
 * nothing but its arguments, so writers can encode on several threads. */
struct LikChunk {
    std::string bytes;
    uint32_t bases = 0, listed = 0, gz = 0, iz = 0, ez = 0, graw = 0, eraw = 0, crc = 0;
    uint8_t cg = 0, ci = 0, ce = 0;
};
LikChunk encodeLikChunk(const LikCode *codes, size_t n, const LikDictionary &dict);

class LikFile {
  public:
    /* Throws std::runtime_error, naming `what`; with ownIo the file closes io. */
    static std::unique_ptr<LikFile> open(lmap_io *io, bool ownIo, const std::string &what);
    ~LikFile();
    LikFile(const LikFile &) = delete;
    LikFile &operator=(const LikFile &) = delete;

    const std::map<std::string, std::string> &meta() const { return _meta; }
    std::string meta(const std::string &key) const;   /* "" if absent */
    int stepsPerDecade() const { return _q; }
    const std::vector<LikCode> &dictionary() const { return _dict; }

    uint32_t numSequences() const { return (uint32_t)_seqs.size(); }
    const std::string &sequenceName(uint32_t i) const { return _seqs[i].name; }
    int64_t sequenceLength(uint32_t i) const { return _seqs[i].length; }

    /* The codes of bases [start, end) of sequence seq.  Thread-safe with a cache per
     * thread (or none).  Throws on a damaged file. */
    void codes(uint32_t seq, int64_t start, int64_t end, std::vector<LikCode> &out, LikCache *cache = nullptr) const;
    /* Their vectors, 4 doubles per base, given the genome's .2bit bases for the same
     * interval (forward); reversed: the reverse strand, last base first, complemented. */
    void probabilities(uint32_t seq, int64_t start, int64_t end, const char *dna, bool reversed, double *out,
                       LikCache *cache = nullptr) const;
    /* Every chunk read, checked and decoded; the directory and dictionary checked.
     * Problems, none if the file is sound. */
    std::vector<std::string> verify() const;

  private:
    LikFile() {}
    struct Seq {
        std::string name;
        int64_t length;
        uint64_t firstChunk;
    };
    lmap_io *_io = nullptr;
    bool _ownIo = false;
    std::string _what;
    int _q = 10;
    uint32_t _chunkBases = 0;
    std::vector<LikCode> _dict;
    std::map<std::string, std::string> _meta;
    std::vector<Seq> _seqs;
    uint64_t _nChunks = 0, _dirOffset = 0;
    uint32_t _dirCrc = 0;
    void decodeChunk(uint64_t chunk, int64_t bases, std::vector<LikCode> &out) const;
};

/* Writing a .hlk: the genome's sequences and the dictionary declared up front, then its
 * chunks in order -- encoded (addChunk), or as codes (add), sequences in order, each
 * whole -- then close().  Writes path.partial and renames it into place on close(); a
 * writer destroyed unclosed removes it. */
class LikWriter {
  public:
    LikWriter(const std::string &path, const std::vector<std::pair<std::string, int64_t>> &seqs,
              const std::map<std::string, std::string> &meta, const LikDictionary &dict, int q = 10,
              uint32_t chunkBases = 65536);
    ~LikWriter();
    LikWriter(const LikWriter &) = delete;
    LikWriter &operator=(const LikWriter &) = delete;
    uint32_t chunkBases() const { return _chunkBases; }
    /* The next chunk: exactly the file's next chunk's bases (chunkBases, or a
     * sequence's rest). */
    void addChunk(const LikChunk &c);
    void add(const LikCode *codes, size_t n);
    void close();

  private:
    std::string _path, _partial;
    FILE *_f = nullptr;
    std::vector<std::pair<std::string, int64_t>> _seqs;
    std::map<std::string, std::string> _meta;
    LikDictionary _dict;
    int _q;
    uint32_t _chunkBases;
    std::vector<LikCode> _chunk;            /* add(): the chunk being filled */
    size_t _seq = 0;                        /* the sequence being written */
    int64_t _inSeq = 0;                     /* its bases so far */
    uint64_t _offset = 64;
    std::string _dir;                       /* directory entries so far */
    uint64_t _nChunks = 0;
    std::vector<uint64_t> _firstChunk;
    int64_t nextChunkBases() const;
    void put(const void *p, size_t n);
};
}

#endif
