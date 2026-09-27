/*
 * HAL format 3: stores -- directory and pack.  See hal3Store.h.
 */
#include "hal3Store.h"
#include "hal3Manifest.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <sstream>
#include <stdexcept>
#include <sys/file.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

using namespace std;

namespace hal3 {

const char PACK_MAGIC[8] = {'\x89', 'H', 'A', 'L', '\r', '\n', '\x1a', '\n'};
static const uint32_t PACK_FORMAT = 3, PACK_VERSION = 1;
static const size_t HEADER_BYTES = 512, SLOT_AT = 16, SLOT_BYTES = 32;
static const char *INDEX_MAGIC = "#hal3-pack-index 1";

namespace {

bool isUrl(const string &p) { return p.find("://") != string::npos; }

lmap_io *openPathIo(const string &path) {
    if (isUrl(path)) {
        lmap_io *u = lmap_io_open_url(path.c_str());
        lmap_io *io = u ? lmap_io_cache(u, 0, 0, 0, 1) : nullptr;
        if (u && !io) lmap_io_close(u);
        return io;
    }
    return lmap_io_open_file(path.c_str());
}

uint64_t get64(const unsigned char *p) {
    uint64_t v = 0;
    for (int k = 7; k >= 0; k--) v = (v << 8) | p[k];
    return v;
}
uint32_t get32(const unsigned char *p) {
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}
void put64(unsigned char *p, uint64_t v) {
    for (int k = 0; k < 8; k++) p[k] = (unsigned char)(v >> (8 * k));
}
void put32(unsigned char *p, uint32_t v) {
    for (int k = 0; k < 4; k++) p[k] = (unsigned char)(v >> (8 * k));
}
uint32_t crcOf(const void *p, size_t n) { return (uint32_t)crc32(crc32(0L, Z_NULL, 0), (const Bytef *)p, (uInt)n); }

string indexText(const map<string, PackEntry> &entries) {
    string s = string(INDEX_MAGIC) + "\n";
    char buf[128];
    for (const auto &e : entries) {
        snprintf(buf, sizeof buf, "\t%llu\t%llu\t%08x\n", (unsigned long long)e.second.offset,
                 (unsigned long long)e.second.length, e.second.crc);
        s += e.first + buf;
    }
    return s;
}

bool parseIndex(const string &text, map<string, PackEntry> &out) {
    istringstream in(text);
    string line;
    if (!getline(in, line) || line != INDEX_MAGIC) return false;
    while (getline(in, line)) {
        if (line.empty()) continue;
        size_t a = line.find('\t'), b = line.find('\t', a + 1), c = line.find('\t', b + 1);
        if (a == string::npos || b == string::npos || c == string::npos) return false;
        PackEntry e;
        char *end;
        e.offset = strtoull(line.c_str() + a + 1, &end, 10);
        e.length = strtoull(line.c_str() + b + 1, &end, 10);
        e.crc = (uint32_t)strtoul(line.c_str() + c + 1, &end, 16);
        out[line.substr(0, a)] = e;
    }
    return true;
}

/* The committed state of the pack on io: the valid slot with the highest generation
 * whose index checks out.  slotOut: which slot that is. */
PackInfo readPack(lmap_io *io, const string &what, int *slotOut) {
    unsigned char hdr[HEADER_BYTES];
    int64_t size = lmap_io_size(io);
    if (size < (int64_t)HEADER_BYTES || lmap_pread(io, hdr, 0, HEADER_BYTES) != 0 || memcmp(hdr, PACK_MAGIC, 8) != 0)
        throw runtime_error(what + " is not a HAL format 3 pack");
    if (get32(hdr + 8) != PACK_FORMAT || get32(hdr + 12) != PACK_VERSION)
        throw runtime_error(what + ": unknown pack version " + to_string(get32(hdr + 8)) + "." + to_string(get32(hdr + 12)));
    int order[2] = {0, 1};
    uint64_t gen[2];
    bool valid[2];
    for (int s = 0; s < 2; s++) {
        const unsigned char *p = hdr + SLOT_AT + SLOT_BYTES * s;
        gen[s] = get64(p);
        valid[s] = gen[s] > 0 && crcOf(p, 28) == get32(p + 28);
    }
    if (gen[1] > gen[0]) swap(order[0], order[1]);
    for (int s : order) {
        if (!valid[s]) continue;
        const unsigned char *p = hdr + SLOT_AT + SLOT_BYTES * s;
        PackInfo info;
        info.generation = gen[s];
        info.indexOffset = get64(p + 8);
        info.indexLength = get64(p + 16);
        info.fileBytes = (uint64_t)size;
        if (info.indexOffset + info.indexLength > (uint64_t)size || info.indexLength > ((uint64_t)1 << 32)) continue;
        string text(info.indexLength, '\0');
        if (lmap_pread(io, &text[0], (int64_t)info.indexOffset, (int64_t)info.indexLength) != 0) continue;
        if (crcOf(text.data(), text.size()) != get32(p + 24) || !parseIndex(text, info.entries)) continue;
        bool inside = true;
        for (auto &e : info.entries) inside = inside && e.second.offset + e.second.length <= info.indexOffset;
        if (!inside) continue;
        if (slotOut) *slotOut = s;
        return info;
    }
    throw runtime_error(what + ": no committed version (an unfinished or damaged pack)");
}

class PackStore : public Store {
  public:
    PackStore(const string &path, lmap_io *io) : _io(io) {
        _path = path;
        _info = readPack(io, path, nullptr);
    }
    ~PackStore() override { lmap_io_close(_io); }
    bool isPack() const override { return true; }
    vector<string> members() const override {
        vector<string> v;
        for (auto &e : _info.entries) v.push_back(e.first);
        return v;
    }
    bool has(const string &name) const override { return _info.entries.count(name) > 0; }
    int64_t size(const string &name) const override { return (int64_t)entry(name).length; }
    lmap_io *openMember(const string &name) const override {
        const PackEntry &e = entry(name);
        lmap_io *io = lmap_io_slice(_io, (int64_t)e.offset, (int64_t)e.length);
        if (!io) throw runtime_error("cannot open " + name + " in " + _path);
        return io;
    }
    const PackInfo &info() const { return _info; }

  private:
    lmap_io *_io;
    PackInfo _info;
    const PackEntry &entry(const string &name) const {
        auto it = _info.entries.find(name);
        if (it == _info.entries.end()) throw runtime_error(_path + " has no member " + name);
        return it->second;
    }
};

class DirStore : public Store {
  public:
    explicit DirStore(const string &dir) {
        _path = dir;
        while (_path.size() > 1 && _path.back() == '/') _path.pop_back();
    }
    bool isPack() const override { return false; }
    vector<string> members() const override {
        vector<string> v;
        if (isUrl(_path)) {
            Manifest m = Manifest::parse(read(MANIFEST_NAME), _path);
            v.push_back(MANIFEST_NAME);
            for (auto &g : m.genomes) {
                v.push_back(g.hgnPath);
                if (!g.dnaPath.empty()) v.push_back(g.dnaPath);
            }
        } else {
            if (has(MANIFEST_NAME)) v.push_back(MANIFEST_NAME);
            if (DIR *d = opendir((_path + "/genomes").c_str())) {
                while (struct dirent *e = readdir(d)) {
                    string n = e->d_name;
                    if (n != "." && n != "..") v.push_back("genomes/" + n);
                }
                closedir(d);
            }
            sort(v.begin(), v.end());
        }
        return v;
    }
    bool has(const string &name) const override {
        if (isUrl(_path)) {
            lmap_io *io = openPathIo(_path + "/" + name);
            if (io) lmap_io_close(io);
            return io != nullptr;
        }
        struct stat st;
        return stat((_path + "/" + name).c_str(), &st) == 0 && S_ISREG(st.st_mode);
    }
    int64_t size(const string &name) const override {
        lmap_io *io = openMember(name);
        int64_t n = lmap_io_size(io);
        lmap_io_close(io);
        return n;
    }
    lmap_io *openMember(const string &name) const override {
        lmap_io *io = openPathIo(_path + "/" + name);
        if (!io) throw runtime_error("cannot open " + _path + "/" + name);
        return io;
    }
};
}

unique_ptr<Store> Store::open(const string &path) {
    if (isUrl(path)) {
        /* a pack, or else a directory */
        if (lmap_io *io = openPathIo(path)) {
            char head[8];
            if (lmap_io_size(io) >= (int64_t)HEADER_BYTES && lmap_pread(io, head, 0, 8) == 0 && !memcmp(head, PACK_MAGIC, 8))
                return unique_ptr<Store>(new PackStore(path, io));
            lmap_io_close(io);
        }
        return unique_ptr<Store>(new DirStore(path));
    }
    struct stat st;
    if (stat(path.c_str(), &st) != 0) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    if (S_ISDIR(st.st_mode)) return unique_ptr<Store>(new DirStore(path));
    lmap_io *io = lmap_io_open_file(path.c_str());
    if (!io) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    char head[16] = {0};
    int64_t n = min<int64_t>(lmap_io_size(io), (int64_t)sizeof head);
    if (n > 0 && lmap_pread(io, head, 0, n) != 0) n = 0;
    if (n >= 8 && !memcmp(head, PACK_MAGIC, 8)) return unique_ptr<Store>(new PackStore(path, io));
    lmap_io_close(io);
    if (Manifest::isManifest(string(head, (size_t)n))) {
        size_t slash = path.rfind('/');
        return unique_ptr<Store>(new DirStore(slash == string::npos ? "." : path.substr(0, slash)));
    }
    throw runtime_error(path + " is not a HAL format 3 alignment (a pack, or a directory with a manifest)");
}

string Store::read(const string &name) const {
    lmap_io *io = openMember(name);
    int64_t n = lmap_io_size(io);
    string s((size_t)max<int64_t>(n, 0), '\0');
    int rc = n > 0 ? lmap_pread(io, &s[0], 0, n) : 0;
    lmap_io_close(io);
    if (rc != 0) throw runtime_error("cannot read " + name + " of " + _path);
    return s;
}

uint32_t Store::crc(const string &name) const {
    lmap_io *io = openMember(name);
    int64_t n = lmap_io_size(io);
    vector<unsigned char> buf((size_t)1 << 20);
    uLong c = crc32(0L, Z_NULL, 0);
    int rc = 0;
    for (int64_t at = 0; at < n && rc == 0; at += (int64_t)buf.size()) {
        int64_t k = min<int64_t>((int64_t)buf.size(), n - at);
        rc = lmap_pread(io, buf.data(), at, k);
        c = crc32(c, buf.data(), (uInt)k);
    }
    lmap_io_close(io);
    if (rc != 0) throw runtime_error("cannot read " + name + " of " + _path);
    return (uint32_t)c;
}

const PackInfo &packInfo(const Store &store) {
    const PackStore *p = dynamic_cast<const PackStore *>(&store);
    if (!p) throw runtime_error(store.path() + " is a directory, not a pack");
    return p->info();
}

/* ---------------------------------------------------------------- writing */

unique_ptr<PackWriter> PackWriter::create(const string &path) {
    unique_ptr<PackWriter> w(new PackWriter());
    w->_path = path;
    w->_writePath = path + ".partial";
    w->_fd = ::open(w->_writePath.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (w->_fd < 0) throw runtime_error("cannot write " + w->_writePath + ": " + strerror(errno));
    if (flock(w->_fd, LOCK_EX | LOCK_NB) != 0) throw runtime_error(w->_writePath + " is being written by another process");
    unsigned char hdr[HEADER_BYTES] = {0};
    memcpy(hdr, PACK_MAGIC, 8);
    put32(hdr + 8, PACK_FORMAT);
    put32(hdr + 12, PACK_VERSION);
    w->append(hdr, sizeof hdr, 0);
    w->_end = HEADER_BYTES;
    w->_renamePending = true;
    return w;
}

unique_ptr<PackWriter> PackWriter::update(const string &path) {
    unique_ptr<PackWriter> w(new PackWriter());
    w->_path = w->_writePath = path;
    w->_fd = ::open(path.c_str(), O_RDWR);
    if (w->_fd < 0) throw runtime_error("cannot open " + path + " for writing: " + strerror(errno));
    if (flock(w->_fd, LOCK_EX | LOCK_NB) != 0) throw runtime_error(path + " is being changed by another process");
    lmap_io *io = lmap_io_open_file(path.c_str());
    if (!io) throw runtime_error("cannot read " + path);
    int slot = 0;
    PackInfo info;
    try {
        info = readPack(io, path, &slot);
    } catch (...) {
        lmap_io_close(io);
        throw;
    }
    lmap_io_close(io);
    w->_entries = w->_opened = info.entries;
    w->_generation = info.generation;
    w->_nextSlot = 1 - slot;
    struct stat st;
    if (fstat(w->_fd, &st) != 0) throw runtime_error("cannot stat " + path);
    w->_end = (uint64_t)st.st_size;   /* after anything an interrupted change left */
    return w;
}

PackWriter::~PackWriter() {
    if (_fd >= 0) ::close(_fd);
    if (_renamePending) ::remove(_writePath.c_str());   /* never committed */
}

void PackWriter::append(const void *data, size_t n, uint64_t at) {
    const char *p = (const char *)data;
    while (n) {
        ssize_t k = pwrite(_fd, p, n, (off_t)at);
        if (k < 0) {
            if (errno == EINTR) continue;
            throw runtime_error("write to " + _writePath + " failed: " + strerror(errno));
        }
        p += k;
        n -= (size_t)k;
        at += (uint64_t)k;
    }
}

uint32_t PackWriter::putIo(const string &name, lmap_io *io) {
    int64_t n = lmap_io_size(io);
    vector<unsigned char> buf((size_t)4 << 20);
    uLong c = crc32(0L, Z_NULL, 0);
    uint64_t start = _end;
    for (int64_t at = 0; at < n; at += (int64_t)buf.size()) {
        int64_t k = min<int64_t>((int64_t)buf.size(), n - at);
        if (lmap_pread(io, buf.data(), at, k) != 0) throw runtime_error("cannot read the source of " + name);
        c = crc32(c, buf.data(), (uInt)k);
        append(buf.data(), (size_t)k, _end);
        _end += (uint64_t)k;
    }
    _entries[name] = PackEntry{start, (uint64_t)n, (uint32_t)c};
    return (uint32_t)c;
}

uint32_t PackWriter::putFile(const string &name, const string &path) {
    lmap_io *io = lmap_io_open_file(path.c_str());
    if (!io) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    uint32_t c;
    try {
        c = putIo(name, io);
    } catch (...) {
        lmap_io_close(io);
        throw;
    }
    lmap_io_close(io);
    return c;
}

uint32_t PackWriter::putBytes(const string &name, const string &bytes) {
    uint64_t start = _end;
    append(bytes.data(), bytes.size(), _end);
    _end += bytes.size();
    uint32_t c = crcOf(bytes.data(), bytes.size());
    _entries[name] = PackEntry{start, (uint64_t)bytes.size(), c};
    return c;
}

void PackWriter::remove(const string &name) { _entries.erase(name); }

void PackWriter::link(const string &from, const string &to) {
    auto it = _opened.find(from);
    if (it == _opened.end()) throw runtime_error(_path + " had no member " + from);
    _entries[to] = it->second;
}

vector<string> PackWriter::members() const {
    vector<string> v;
    for (auto &e : _entries) v.push_back(e.first);
    return v;
}

void PackWriter::commit() {
    string idx = indexText(_entries);
    uint64_t at = _end;
    append(idx.data(), idx.size(), at);
    _end += idx.size();
    if (fsync(_fd) != 0) throw runtime_error("fsync of " + _writePath + " failed: " + strerror(errno));
    unsigned char slot[SLOT_BYTES];
    put64(slot, ++_generation);
    put64(slot + 8, at);
    put64(slot + 16, idx.size());
    put32(slot + 24, crcOf(idx.data(), idx.size()));
    put32(slot + 28, crcOf(slot, 28));
    append(slot, sizeof slot, SLOT_AT + SLOT_BYTES * (uint64_t)_nextSlot);
    if (fsync(_fd) != 0) throw runtime_error("fsync of " + _writePath + " failed: " + strerror(errno));
    _nextSlot = 1 - _nextSlot;
    if (_renamePending) {
        if (::rename(_writePath.c_str(), _path.c_str()) != 0)
            throw runtime_error("cannot rename " + _writePath + " to " + _path + ": " + strerror(errno));
        _renamePending = false;
        _writePath = _path;
    }
}

uint64_t packAlignment(const string &in, const string &out) {
    unique_ptr<Store> store = Store::open(in);
    Manifest m = Manifest::parse(store->read(MANIFEST_NAME), in);
    vector<pair<string, const uint32_t *>> named;
    for (const GenomeEntry &g : m.genomes) {
        named.push_back({g.hgnPath, &g.hgnCrc});
        if (!g.dnaPath.empty()) named.push_back({g.dnaPath, &g.dnaCrc});
    }
    named.push_back({MANIFEST_NAME, nullptr});
    unique_ptr<PackWriter> w = PackWriter::create(out);
    uint64_t bytes = 0;
    for (auto &n : named) {
        lmap_io *io = store->openMember(n.first);
        uint32_t c;
        try {
            c = w->putIo(n.first, io);
        } catch (...) {
            lmap_io_close(io);
            throw;
        }
        bytes += (uint64_t)lmap_io_size(io);
        lmap_io_close(io);
        if (n.second && c != *n.second) throw runtime_error(n.first + " does not match the manifest's CRC-32; nothing written");
    }
    w->commit();
    return bytes;
}
}
