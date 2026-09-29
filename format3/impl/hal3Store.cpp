/*
 * HAL format 3: stores -- directory and pack.  See hal3Store.h.
 */
#include "hal3Store.h"
#include "hal3Manifest.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <deque>
#include <list>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <sys/file.h>
#include <sys/resource.h>
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

/* Descriptors for a local directory's members, at most `cap` open at a time (see
 * hal3Store.h).  A read takes its member's descriptor -- opening it, after closing the
 * least recently used idle one if the pool is full -- and gives it back after.  If
 * every descriptor is busy, a read opens one of its own and closes it after.
 * Thread-safe: the bookkeeping is under a lock; the reads are not. */
class FilePool {
  public:
    explicit FilePool(size_t cap) : _cap(max<size_t>(1, cap)) {}
    ~FilePool() {
        for (File &f : _files)
            if (f.fd >= 0) ::close(f.fd);
    }
    FilePool(const FilePool &) = delete;
    FilePool &operator=(const FilePool &) = delete;

    /* The member at path (one id per path, however often it is opened). */
    size_t add(const string &path) {
        lock_guard<mutex> lock(_m);
        auto i = _ids.find(path);
        if (i != _ids.end()) return i->second;
        _files.emplace_back();
        _files.back().path = path;
        return _ids[path] = _files.size() - 1;
    }
    int64_t read(size_t id, void *buf, int64_t off, int64_t len) {
        int fd = acquire(id);
        if (fd < 0) return -1;
        int64_t got = 0;
        bool failed = false;
        while (got < len) {
            ssize_t k = ::pread(fd, (char *)buf + got, (size_t)(len - got), (off_t)(off + got));
            if (k < 0 && errno == EINTR) continue;
            if (k < 0) failed = true;
            if (k <= 0) break;   /* an error, or the end of the file */
            got += k;
        }
        release(id);
        return failed ? -1 : got;
    }

  private:
    struct File {
        string path;
        int fd = -1;
        int users = 0;
        bool idle = false;                 /* in _idle */
        list<size_t>::iterator at;
    };
    mutex _m;
    const size_t _cap;
    size_t _open = 0;
    deque<File> _files;
    map<string, size_t> _ids;
    list<size_t> _idle;                   /* open and unused, most recently used first */

    void closeIdle() {                    /* with _m held */
        size_t v = _idle.back();
        _idle.pop_back();
        File &x = _files[v];
        x.idle = false;
        ::close(x.fd);
        x.fd = -1;
        _open--;
    }
    int acquire(size_t id) {
        lock_guard<mutex> lock(_m);
        File &f = _files[id];
        if (f.fd >= 0) {
            if (f.users++ == 0 && f.idle) {
                _idle.erase(f.at);
                f.idle = false;
            }
            return f.fd;
        }
        while (_open >= _cap && !_idle.empty()) closeIdle();
        int fd = ::open(f.path.c_str(), O_RDONLY | O_CLOEXEC);
        while (fd < 0 && (errno == EMFILE || errno == ENFILE) && !_idle.empty()) {   /* others' files count too */
            closeIdle();
            fd = ::open(f.path.c_str(), O_RDONLY | O_CLOEXEC);
        }
        if (fd < 0) return -1;
        f.fd = fd;
        f.users = 1;
        _open++;
        return fd;
    }
    void release(size_t id) {
        lock_guard<mutex> lock(_m);
        File &f = _files[id];
        if (--f.users > 0) return;
        if (_open > _cap) {               /* opened past the cap while all were busy */
            ::close(f.fd);
            f.fd = -1;
            _open--;
            return;
        }
        _idle.push_front(id);
        f.at = _idle.begin();
        f.idle = true;
    }
};

/* A member of a local directory as an lmap_io over the pool. */
struct PooledMember {
    shared_ptr<FilePool> pool;
    size_t id;
};
int64_t pooledRead(void *ctx, void *buf, int64_t off, int64_t len) {
    PooledMember *m = (PooledMember *)ctx;
    return m->pool->read(m->id, buf, off, len);
}
void pooledClose(void *ctx) { delete (PooledMember *)ctx; }
const lmap_io_backend POOLED = {"hal3-pool", pooledRead, pooledClose};

size_t defaultMaxOpenFiles() {
    if (const char *e = getenv("HAL3_MAX_OPEN_FILES"))
        if (atoll(e) > 0) return (size_t)atoll(e);
    struct rlimit rl;
    size_t soft = getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY ? (size_t)rl.rlim_cur : 4096;
    return max<size_t>(4, min<size_t>(256, soft / 4));
}

class DirStore : public Store {
  public:
    DirStore(const string &dir, size_t maxOpenFiles)
        : _pool(new FilePool(maxOpenFiles ? maxOpenFiles : defaultMaxOpenFiles())) {
        _path = dir;
        while (_path.size() > 1 && _path.back() == '/') _path.pop_back();
    }
    bool isPack() const override { return false; }
    vector<string> members() const override {
        vector<string> v;
        if (isUrl(_path)) {
            Manifest m = Manifest::parse(read(MANIFEST_NAME), _path);
            v.push_back(MANIFEST_NAME);
            for (const NamedFile &f : m.files()) v.push_back(f.path);
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
        const string p = _path + "/" + name;
        if (isUrl(_path)) {
            lmap_io *io = openPathIo(p);
            if (!io) throw runtime_error("cannot open " + p);
            return io;
        }
        struct stat st;
        if (stat(p.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) throw runtime_error("cannot open " + p + ": " + strerror(errno));
        int fd = ::open(p.c_str(), O_RDONLY | O_CLOEXEC);   /* readable now, as a plain open would say */
        if (fd < 0) throw runtime_error("cannot open " + p + ": " + strerror(errno));
        ::close(fd);
        PooledMember *m = new PooledMember{_pool, _pool->add(p)};
        lmap_io *io = lmap_io_open_backend(&POOLED, m, (int64_t)st.st_size);
        if (!io) {
            delete m;
            throw runtime_error("cannot open " + p);
        }
        return io;
    }

  private:
    shared_ptr<FilePool> _pool;
};
}

unique_ptr<Store> Store::open(const string &path, size_t maxOpenFiles) {
    if (isUrl(path)) {
        /* a pack, or else a directory */
        if (lmap_io *io = openPathIo(path)) {
            char head[8];
            if (lmap_io_size(io) >= (int64_t)HEADER_BYTES && lmap_pread(io, head, 0, 8) == 0 && !memcmp(head, PACK_MAGIC, 8))
                return unique_ptr<Store>(new PackStore(path, io));
            lmap_io_close(io);
        }
        return unique_ptr<Store>(new DirStore(path, maxOpenFiles));
    }
    struct stat st;
    if (stat(path.c_str(), &st) != 0) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    if (S_ISDIR(st.st_mode)) return unique_ptr<Store>(new DirStore(path, maxOpenFiles));
    lmap_io *io = lmap_io_open_file(path.c_str());
    if (!io) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    char head[16] = {0};
    int64_t n = min<int64_t>(lmap_io_size(io), (int64_t)sizeof head);
    if (n > 0 && lmap_pread(io, head, 0, n) != 0) n = 0;
    if (n >= 8 && !memcmp(head, PACK_MAGIC, 8)) return unique_ptr<Store>(new PackStore(path, io));
    lmap_io_close(io);
    if (Manifest::isManifest(string(head, (size_t)n))) {
        size_t slash = path.rfind('/');
        return unique_ptr<Store>(new DirStore(slash == string::npos ? "." : path.substr(0, slash), maxOpenFiles));
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
    vector<NamedFile> files = m.files();
    vector<pair<string, const uint32_t *>> named;
    for (const NamedFile &f : files) named.push_back({f.path, &f.crc});
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

/* Copy an lmap_io to a file; returns the CRC-32 and sets *bytes. */
uint32_t copyToFile(lmap_io *io, const string &path, uint64_t *bytes) {
    FILE *f = fopen(path.c_str(), "wb");
    if (!f) throw runtime_error("cannot write " + path + ": " + strerror(errno));
    int64_t n = lmap_io_size(io);
    vector<unsigned char> buf((size_t)4 << 20);
    uLong c = crc32(0L, Z_NULL, 0);
    bool good = true;
    for (int64_t at = 0; at < n && good; at += (int64_t)buf.size()) {
        int64_t k = min<int64_t>((int64_t)buf.size(), n - at);
        good = lmap_pread(io, buf.data(), at, k) == 0 && fwrite(buf.data(), 1, (size_t)k, f) == (size_t)k;
        c = crc32(c, buf.data(), (uInt)k);
    }
    good = (fclose(f) == 0) && good;
    if (!good) throw runtime_error("cannot write " + path);
    *bytes = (uint64_t)n;
    return (uint32_t)c;
}
}
