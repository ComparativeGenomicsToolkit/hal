/*
 * hal pack / unpack / ls: between the single-file (pack) and directory forms of an
 * alignment, and what a store holds.
 */
#include "hal3Cli.h"
#include "hal3Manifest.h"
#include "hal3Store.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <sys/stat.h>
#include <vector>

using namespace std;

namespace hal3 {

namespace {

/* Options: only --help. */
void noOptions(int argc, char **argv, const char *usage) {
    static const option longOpts[] = {{"help", no_argument, nullptr, 'h'}, {nullptr, 0, nullptr, 0}};
    int c;
    while ((c = getopt_long(argc, argv, "h", longOpts, nullptr)) != -1) {
        if (c == 'h') {
            cout << usage;
            exit(0);
        }
        throw CliError("bad option");
    }
}

/* The members an alignment is made of, as its manifest names them (manifest last), with
 * the CRC-32 the manifest records for each (the manifest's own: none). */
vector<pair<string, const uint32_t *>> named(const Manifest &m) {
    vector<pair<string, const uint32_t *>> v;
    for (const GenomeEntry &g : m.genomes) {
        v.push_back({g.hgnPath, &g.hgnCrc});
        if (!g.dnaPath.empty()) v.push_back({g.dnaPath, &g.dnaCrc});
    }
    v.push_back({MANIFEST_NAME, nullptr});
    return v;
}
}

static const char *PACK_USAGE =
    "usage: hal pack <alignment> <out.hal>\n"
    "Write an alignment -- a directory, or a pack -- as a single file (a pack).  Only what\n"
    "the manifest names is copied, each member checked against the manifest's CRC-32, so\n"
    "packing a pack drops the dead space its updates left behind.  out.hal may be the input:\n"
    "it is written beside it and renamed into place.\n";

int main_pack(int argc, char **argv) {
    return runCommand("pack", [&]() {
        noOptions(argc, argv, PACK_USAGE);
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        string in = argv[optind], out = argv[optind + 1];
        uint64_t bytes = packAlignment(in, out);
        Manifest m = Manifest::parse(Store::open(out)->read(MANIFEST_NAME), out);
        cerr << "hal pack: " << out << ": " << m.genomes.size() << " genomes, " << bytes << " bytes" << endl;
        return 0;
    });
}

static const char *UNPACK_USAGE =
    "usage: hal unpack <alignment> <outDir>\n"
    "Write an alignment as a directory: outDir/manifest and outDir/genomes/.  outDir must not\n"
    "already hold an alignment.\n";

int main_unpack(int argc, char **argv) {
    return runCommand("unpack", [&]() {
        noOptions(argc, argv, UNPACK_USAGE);
        if (argc - optind != 2) throw CliError("expected 2 arguments, got " + to_string(argc - optind));
        string in = argv[optind], dir = argv[optind + 1];
        unique_ptr<Store> store = Store::open(in);
        Manifest m = Manifest::parse(store->read(MANIFEST_NAME), in);
        struct stat st;
        if (stat((dir + "/" + MANIFEST_NAME).c_str(), &st) == 0) throw runtime_error(dir + " already holds an alignment");
        for (const string &d : {dir, dir + "/genomes"})
            if (mkdir(d.c_str(), 0777) != 0 && errno != EEXIST) throw runtime_error("cannot create " + d + ": " + strerror(errno));
        vector<char> buf((size_t)4 << 20);
        for (auto &n : named(m)) {
            /* the manifest last, so an interrupted unpack is not taken for a finished one */
            string path = dir + "/" + n.first, partial = path + ".partial";
            lmap_io *io = store->openMember(n.first);
            FILE *f = fopen(partial.c_str(), "wb");
            if (!f) {
                lmap_io_close(io);
                throw runtime_error("cannot write " + partial + ": " + strerror(errno));
            }
            int64_t size = lmap_io_size(io);
            bool ok = true;
            for (int64_t at = 0; at < size && ok; at += (int64_t)buf.size()) {
                int64_t k = min<int64_t>((int64_t)buf.size(), size - at);
                ok = lmap_pread(io, buf.data(), at, k) == 0 && fwrite(buf.data(), 1, (size_t)k, f) == (size_t)k;
            }
            lmap_io_close(io);
            ok = (fclose(f) == 0) && ok;
            if (!ok || rename(partial.c_str(), path.c_str()) != 0) {
                remove(partial.c_str());
                throw runtime_error("cannot write " + path);
            }
        }
        cerr << "hal unpack: " << dir << ": " << m.genomes.size() << " genomes" << endl;
        return 0;
    });
}

static const char *LS_USAGE =
    "usage: hal ls <alignment>\n"
    "The members of an alignment and their sizes; for a pack, where each member lies, its\n"
    "CRC-32, and how much of the file is dead space left by updates.\n";

int main_ls(int argc, char **argv) {
    return runCommand("ls", [&]() {
        noOptions(argc, argv, LS_USAGE);
        if (argc - optind != 1) throw CliError("expected 1 argument, got " + to_string(argc - optind));
        unique_ptr<Store> store = Store::open(argv[optind]);
        if (store->isPack()) {
            const PackInfo &info = packInfo(*store);
            uint64_t live = 512 + info.indexLength;
            for (auto &e : info.entries) live += e.second.length;
            printf("# pack %s: generation %llu, %llu bytes, %llu live, %llu dead (%.1f%%)\n", store->path().c_str(),
                   (unsigned long long)info.generation, (unsigned long long)info.fileBytes, (unsigned long long)live,
                   (unsigned long long)(info.fileBytes - live), info.fileBytes ? 100.0 * (info.fileBytes - live) / info.fileBytes : 0.0);
            printf("#member\toffset\tlength\tcrc32\n");
            for (auto &e : info.entries)
                printf("%s\t%llu\t%llu\t%08x\n", e.first.c_str(), (unsigned long long)e.second.offset,
                       (unsigned long long)e.second.length, e.second.crc);
        } else {
            printf("# directory %s\n#member\tlength\n", store->path().c_str());
            for (const string &n : store->members()) printf("%s\t%lld\n", n.c_str(), (long long)store->size(n));
        }
        if (fflush(stdout) != 0) throw runtime_error("write failed");
        return 0;
    });
}
}
