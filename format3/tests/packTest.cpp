/*
 * packTest: the pack's commit protocol under simulated crashes and damage.  Takes a copy
 * of a committed pack and changes it: an update abandoned before commit, a damaged commit
 * slot, a torn index, an interrupted create, and two writers at once.
 *
 * usage: packTest pack.hal   (the file is modified)
 */
#include "hal3Store.h"
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
using namespace std;
using namespace hal3;
static int fails = 0;
static void check(const char *what, bool ok) { printf("  %s %s\n", ok ? "ok  " : "FAIL", what); if (!ok) fails++; }
static uint64_t gen(const string &p) { return packInfo(*Store::open(p)).generation; }
static bool hasMember(const string &p, const string &m) { return Store::open(p)->has(m); }
static off_t fsize(const string &p) { struct stat st; stat(p.c_str(), &st); return st.st_size; }
int main(int argc, char **argv) {
    string p = argv[1];                     /* a copy of a pack, generation 1 */
    uint64_t g0 = gen(p);
    string big(10 << 20, 'x');
    {   /* 1: an update abandoned before commit */
        auto w = PackWriter::update(p);
        w->putBytes("genomes/junk", big);
        w->remove("manifest");
    }
    check("abandoned update: same generation, members unchanged", gen(p) == g0 && !hasMember(p, "genomes/junk") && hasMember(p, "manifest"));
    check("abandoned update left its bytes as dead space", fsize(p) > (off_t)(10 << 20));
    {   /* 2: a committed update after that */
        auto w = PackWriter::update(p);
        w->putBytes("genomes/extra", "hello");
        w->commit();
    }
    check("committed update: generation + 1, member present", gen(p) == g0 + 1 && hasMember(p, "genomes/extra"));
    {   /* 3: damage the newest slot: the previous version is current again */
        int fd = open(p.c_str(), O_RDWR);
        unsigned char b[64];
        pread(fd, b, 64, 16);
        uint64_t g1 = 0, g2 = 0;
        for (int k = 7; k >= 0; k--) { g1 = (g1 << 8) | b[k]; g2 = (g2 << 8) | b[32 + k]; }
        int newest = g2 > g1 ? 1 : 0;
        unsigned char x = b[32 * newest + 3] ^ 0xff;
        pwrite(fd, &x, 1, 16 + 32 * newest + 3);
        close(fd);
    }
    check("damaged newest slot: falls back to the previous version", gen(p) == g0 && !hasMember(p, "genomes/extra"));
    {   /* 4: the writer after that commits over the damaged slot's generation */
        auto w = PackWriter::update(p);
        w->putBytes("genomes/extra2", "again");
        w->commit();
    }
    check("update after damage: commits, previous members intact", hasMember(p, "genomes/extra2") && hasMember(p, "manifest"));
    {   /* 5: a torn index: truncate the newest commit's index */
        uint64_t g = gen(p);
        auto w = PackWriter::update(p);
        w->putBytes("genomes/extra3", "torn");
        w->commit();
        w.reset();
        truncate(p.c_str(), fsize(p) - 10);
        check("torn index: falls back to the previous version", gen(p) == g && !hasMember(p, "genomes/extra3"));
    }
    {   /* 6: an interrupted create leaves nothing behind */
        string q = p + ".new";
        {
            auto w = PackWriter::create(q);
            w->putBytes("manifest", "x");
        }
        struct stat st;
        check("interrupted create: no file, no .partial", stat(q.c_str(), &st) != 0 && stat((q + ".partial").c_str(), &st) != 0);
    }
    {   /* 7: two writers at once */
        auto w = PackWriter::update(p);
        bool refused = false;
        try { PackWriter::update(p); } catch (exception &) { refused = true; }
        check("a second writer is refused while one holds the pack", refused);
    }
    printf("%s (%d failures)\n", fails ? "FAILURES" : "ALL PASS", fails);
    return fails ? 1 : 0;
}
