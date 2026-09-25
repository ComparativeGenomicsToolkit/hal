/*
 * HAL format 3: the manifest.  See hal3Manifest.h.
 */
#include "hal3Manifest.h"

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <zlib.h>

using namespace std;

namespace hal3 {

const char *const MANIFEST_MAGIC = "#hal-manifest 3";
const char *const MANIFEST_NAME = "manifest";
const char *const HGN_PROFILE = "hal.branch";

static string escape(const string &s) {
    string o;
    for (char c : s) {
        if (c == '\t') o += "\\t";
        else if (c == '\n') o += "\\n";
        else if (c == '\\') o += "\\\\";
        else o += c;
    }
    return o;
}

static string unescape(const string &s, size_t lineNo) {
    string o;
    for (size_t i = 0; i < s.size(); i++) {
        if (s[i] != '\\') { o += s[i]; continue; }
        if (++i == s.size()) throw runtime_error("manifest line " + to_string(lineNo) + ": trailing backslash");
        if (s[i] == 't') o += '\t';
        else if (s[i] == 'n') o += '\n';
        else if (s[i] == '\\') o += '\\';
        else throw runtime_error("manifest line " + to_string(lineNo) + ": unknown escape \\" + s[i]);
    }
    return o;
}

static vector<string> splitTabs(const string &line) {
    vector<string> f;
    size_t start = 0;
    for (;;) {
        size_t t = line.find('\t', start);
        f.push_back(line.substr(start, t == string::npos ? string::npos : t - start));
        if (t == string::npos) break;
        start = t + 1;
    }
    return f;
}

/* The shortest decimal that reads back as exactly d. */
static string fmtDouble(double d) {
    char buf[64];
    for (int prec = 6; prec <= 17; prec++) {
        snprintf(buf, sizeof buf, "%.*g", prec, d);
        if (strtod(buf, NULL) == d) break;
    }
    return buf;
}

const GenomeEntry *Manifest::find(const string &name) const {
    for (const GenomeEntry &g : genomes)
        if (g.name == name) return &g;
    return NULL;
}

string Manifest::text() const {
    string out = string(MANIFEST_MAGIC) + "\n";
    out += "tree\t" + escape(newick) + "\n";
    for (const auto &kv : meta) out += "meta\t" + escape(kv.first) + "\t" + escape(kv.second) + "\n";
    char buf[256];
    for (const GenomeEntry &g : genomes) {
        out += "genome\t" + escape(g.name) + "\t" + (g.parent.empty() ? string("-") : escape(g.parent)) + "\t" +
               (g.parent.empty() ? string("-") : fmtDouble(g.branchLength));
        snprintf(buf, sizeof buf, "\t%llu\t%llu\t", (unsigned long long)g.length, (unsigned long long)g.numSequences);
        out += buf + escape(g.hgnPath);
        snprintf(buf, sizeof buf, "\t%llu\t%08x\t", (unsigned long long)g.hgnBytes, g.hgnCrc);
        out += buf + (g.dnaPath.empty() ? string("-") : escape(g.dnaPath));
        snprintf(buf, sizeof buf, "\t%llu\t%08x\n", (unsigned long long)g.dnaBytes, g.dnaCrc);
        out += buf;
    }
    return out;
}

void Manifest::write(const string &path) const {
    string partial = path + ".partial", t = text();
    FILE *f = fopen(partial.c_str(), "w");
    if (!f) throw runtime_error("cannot write " + partial + ": " + strerror(errno));
    fwrite(t.data(), 1, t.size(), f);
    if (ferror(f) | fclose(f)) {
        remove(partial.c_str());
        throw runtime_error("cannot write " + partial);
    }
    if (rename(partial.c_str(), path.c_str()) != 0) {
        remove(partial.c_str());
        throw runtime_error("cannot rename " + partial + " to " + path + ": " + strerror(errno));
    }
}

static uint64_t parseU64(const string &s, size_t lineNo, const char *what) {
    char *end = NULL;
    errno = 0;
    unsigned long long v = strtoull(s.c_str(), &end, 10);
    if (s.empty() || *end || errno || s[0] == '-')
        throw runtime_error("manifest line " + to_string(lineNo) + ": bad " + what + " '" + s + "'");
    return v;
}

static uint32_t parseCrc(const string &s, size_t lineNo) {
    char *end = NULL;
    unsigned long v = strtoul(s.c_str(), &end, 16);
    if (s.size() != 8 || *end) throw runtime_error("manifest line " + to_string(lineNo) + ": bad crc32 '" + s + "'");
    return (uint32_t)v;
}

Manifest Manifest::read(const string &path) {
    ifstream f(path);
    if (!f) throw runtime_error("cannot open manifest " + path);
    stringstream ss;
    ss << f.rdbuf();
    return parse(ss.str(), path);
}

Manifest Manifest::parse(const string &text, const string &path) {
    istringstream in(text);
    Manifest m;
    string line;
    size_t lineNo = 0;
    bool haveTree = false;
    while (getline(in, line)) {
        lineNo++;
        if (lineNo == 1) {
            if (line != MANIFEST_MAGIC) throw runtime_error(path + " is not a HAL format 3 manifest");
            continue;
        }
        if (line.empty() || line[0] == '#') continue;
        vector<string> f = splitTabs(line);
        if (f[0] == "tree" && f.size() == 2) {
            if (haveTree) throw runtime_error("manifest line " + to_string(lineNo) + ": second tree");
            m.newick = unescape(f[1], lineNo);
            haveTree = true;
        } else if (f[0] == "meta" && f.size() == 3) {
            m.meta[unescape(f[1], lineNo)] = unescape(f[2], lineNo);
        } else if (f[0] == "genome" && f.size() == 12) {
            GenomeEntry g;
            g.name = unescape(f[1], lineNo);
            if (f[2] != "-") {
                g.parent = unescape(f[2], lineNo);
                char *end = NULL;
                g.branchLength = strtod(f[3].c_str(), &end);
                if (f[3].empty() || *end) throw runtime_error("manifest line " + to_string(lineNo) + ": bad branch length");
            }
            g.length = parseU64(f[4], lineNo, "length");
            g.numSequences = parseU64(f[5], lineNo, "sequence count");
            g.hgnPath = unescape(f[6], lineNo);
            g.hgnBytes = parseU64(f[7], lineNo, "size");
            g.hgnCrc = parseCrc(f[8], lineNo);
            if (f[9] != "-") g.dnaPath = unescape(f[9], lineNo);
            g.dnaBytes = parseU64(f[10], lineNo, "size");
            g.dnaCrc = parseCrc(f[11], lineNo);
            if (m.find(g.name)) throw runtime_error("manifest line " + to_string(lineNo) + ": genome " + g.name + " listed twice");
            m.genomes.push_back(g);
        } else {
            throw runtime_error("manifest line " + to_string(lineNo) + ": unrecognised record '" + f[0] + "'");
        }
    }
    if (lineNo == 0) throw runtime_error(path + " is empty");
    if (!haveTree) throw runtime_error(path + ": no tree record");
    return m;
}

bool Manifest::isManifest(const string &initialBytes) {
    size_t n = strlen(MANIFEST_MAGIC);
    return initialBytes.size() > n && initialBytes.compare(0, n, MANIFEST_MAGIC) == 0 &&
           (initialBytes[n] == '\n' || initialBytes[n] == '\r');
}

string genomeFileStem(const string &name) {
    string o;
    for (size_t i = 0; i < name.size(); i++) {
        unsigned char c = name[i];
        bool ok = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
                  c == '_' || c == '-' || (c == '.' && i > 0);
        if (ok) o += (char)c;
        else {
            char buf[4];
            snprintf(buf, sizeof buf, "%%%02X", c);
            o += buf;
        }
    }
    return o;
}

uint32_t fileCrc32(const string &path, uint64_t *bytes) {
    FILE *f = fopen(path.c_str(), "rb");
    if (!f) throw runtime_error("cannot open " + path + ": " + strerror(errno));
    vector<unsigned char> buf(1 << 20);
    uLong crc = crc32(0L, Z_NULL, 0);
    uint64_t total = 0;
    size_t n;
    while ((n = fread(buf.data(), 1, buf.size(), f)) > 0) {
        crc = crc32(crc, buf.data(), (uInt)n);
        total += n;
    }
    bool err = ferror(f);
    fclose(f);
    if (err) throw runtime_error("read error on " + path);
    if (bytes) *bytes = total;
    return (uint32_t)crc;
}

/* Newick: (child,child,...)name:length ... ;  Names may be quoted with '. */
namespace {
struct NewickParser {
    const string &s;
    size_t p = 0;
    vector<NewickNode> nodes;
    explicit NewickParser(const string &str) : s(str) {}
    [[noreturn]] void fail(const string &why) {
        throw runtime_error("bad Newick tree at offset " + to_string(p) + ": " + why);
    }
    void ws() { while (p < s.size() && isspace((unsigned char)s[p])) p++; }
    string label() {
        ws();
        string o;
        if (p < s.size() && s[p] == '\'') {
            for (p++; p < s.size(); p++) {
                if (s[p] == '\'') {
                    if (p + 1 < s.size() && s[p + 1] == '\'') { o += '\''; p++; }
                    else { p++; return o; }
                } else o += s[p];
            }
            fail("unterminated quoted name");
        }
        while (p < s.size() && !strchr("(),:;", s[p]) && !isspace((unsigned char)s[p])) o += s[p++];
        return o;
    }
    int node(int parent) {
        int me = (int)nodes.size();
        nodes.push_back({"", parent, -1});
        ws();
        if (p < s.size() && s[p] == '(') {
            p++;
            for (;;) {
                node(me);
                ws();
                if (p < s.size() && s[p] == ',') { p++; continue; }
                if (p < s.size() && s[p] == ')') { p++; break; }
                fail("expected ',' or ')'");
            }
        }
        nodes[me].name = label();
        ws();
        if (p < s.size() && s[p] == ':') {
            p++;
            ws();
            const char *b = s.c_str() + p;
            char *end = NULL;
            double v = strtod(b, &end);
            if (end == b) fail("bad branch length");
            p += (size_t)(end - b);
            nodes[me].branchLength = v;
        }
        return me;
    }
};
}

vector<NewickNode> parseNewick(const string &newick) {
    NewickParser np(newick);
    np.node(-1);
    np.ws();
    if (np.p >= newick.size() || newick[np.p] != ';') np.fail("expected ';'");
    np.p++;
    np.ws();
    if (np.p != newick.size()) np.fail("text after ';'");
    return np.nodes;
}
}
