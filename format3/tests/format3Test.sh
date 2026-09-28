#!/bin/bash
# HAL format 3 tests: generated alignments only (halRandGen), so they need nothing
# outside this repository and its sibling liftmap / twobit64 checkouts.
#
#   - halRandGen alignments (hdf5 and mmap back ends, several seeds) convert, pass
#     halToFormat3's own check against the input, and pass hal3Validate
#   - the same alignment gives byte-identical genome files from either back end
#   - an existing alignment is not overwritten
#   - hal validate catches corruption, missing and orphaned files, a swapped genome
#     file, a wrong parent, and a mis-ordered roster
#   - hal convert --checkOnly catches a genome's DNA replaced by another's
#
# usage: format3Test.sh BINDIR
set -u
BIN=${1:-../bin}
TMP=$(mktemp -d -t format3Test.XXXXXX)
trap 'rm -rf "$TMP"' EXIT
fails=0
check() {   # check "what" command...
    local what=$1; shift
    if "$@" >"$TMP/out" 2>&1; then echo "  ok   $what"; else echo "  FAIL $what"; sed 's/^/       /' "$TMP/out" | head -5; fails=$((fails + 1)); fi
}
refuses() {   # refuses "what" command...  (must exit non-zero, and not by crashing)
    local what=$1; shift
    "$@" >"$TMP/out" 2>&1
    local status=$?
    if [ $status -eq 0 ]; then echo "  FAIL $what (was accepted)"; fails=$((fails + 1))
    elif [ $status -ge 128 ]; then echo "  FAIL $what (crashed: exit status $status)"; fails=$((fails + 1))
    else echo "  ok   $what"; fi
}

echo "conversion:"
for seed in 1 7 23; do
    for fmt in hdf5 mmap; do
        "$BIN/halRandGen" --testRand --seed $seed --format $fmt "$TMP/r$seed.$fmt.hal" >/dev/null 2>&1
        check "seed $seed, $fmt: converts and passes the check against HAL" \
            "$BIN/hal" convert "$TMP/r$seed.$fmt.hal" "$TMP/r$seed.$fmt" --threads 3 --check 2000 --directory
        check "seed $seed, $fmt: hal3Validate passes" "$BIN/hal" validate "$TMP/r$seed.$fmt"
    done
    check "seed $seed: hdf5 and mmap inputs give identical genome files" \
        diff -r "$TMP/r$seed.hdf5/genomes" "$TMP/r$seed.mmap/genomes"
done
refuses "an existing alignment is not overwritten" "$BIN/hal" convert "$TMP/r1.hdf5.hal" "$TMP/r1.hdf5" --directory
check "--checkOnly passes on an intact alignment" "$BIN/hal" convert "$TMP/r1.hdf5.hal" "$TMP/r1.hdf5" --checkOnly

echo "hal validate catches:"
A="$TMP/r7.hdf5"
G=($(ls "$A/genomes" | sed -n 's/\.hgn$//p'))       # genome file stems
fresh() { rm -rf "$TMP/t"; cp -r "$A" "$TMP/t"; }
fresh; printf '\x55' | dd of="$TMP/t/genomes/${G[1]}.hgn" bs=1 seek=300 conv=notrunc 2>/dev/null
refuses "a corrupted .hgn" "$BIN/hal" validate "$TMP/t"
fresh; printf '\x55' | dd of="$TMP/t/genomes/${G[1]}.hgn" bs=1 seek=300 conv=notrunc 2>/dev/null
refuses "a corrupted .hgn, even with --noCrc" "$BIN/hal" validate "$TMP/t" --noCrc
fresh; rm "$TMP/t/genomes/${G[2]}.2bit"
refuses "a missing .2bit" "$BIN/hal" validate "$TMP/t"
fresh; cp "$TMP/t/genomes/${G[1]}.hgn" "$TMP/t/genomes/extra.hgn"
refuses "an orphaned genome file" "$BIN/hal" validate "$TMP/t"
fresh; cp "$TMP/t/genomes/${G[1]}.hgn" "$TMP/t/genomes/${G[2]}.hgn"
refuses "one genome's .hgn in place of another's (--noCrc)" "$BIN/hal" validate "$TMP/t" --noCrc
# re-parent the first non-root genome (a child of the root) under the last one listed
fresh; awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {if ($1=="genome") last=$2; next}
    $1=="genome" && $3!="-" && !done {$3=last; done=1} {print}' "$A/manifest" "$A/manifest" >"$TMP/t/manifest"
refuses "a wrong parent in the roster" "$BIN/hal" validate "$TMP/t"
fresh; awk 'NR==1{print; next} /^genome\t/ && !held {held=$0; next} {print} END{print held}' "$A/manifest" >"$TMP/t/manifest"
refuses "the root not listed first" "$BIN/hal" validate "$TMP/t"

echo "hal convert --checkOnly catches:"
fresh; cp "$TMP/t/genomes/${G[1]}.2bit" "$TMP/t/genomes/${G[2]}.2bit"
refuses "one genome's DNA in place of another's" "$BIN/hal" convert "$TMP/r7.hdf5.hal" "$TMP/t" --checkOnly

echo "hal liftover against halLiftover:"
if [ ! -x "$BIN/halLiftover" ] || [ ! -x "$BIN/halStats" ]; then
    echo "  skip (needs halLiftover and halStats: build the liftover and stats modules)"
else
    # With paralogs traversed, the output must be byte-identical to halLiftover's, or --
    # where its line grouping depends on HAL's segmentation, in overlapping paralogs --
    # the same aligned bases, strands and fields.  --noDupes follows HAL's stored
    # canonical-paralog pointers, which format 3 replaces with a rule (the copy with the
    # smallest child coordinate); on Cactus alignments they agree, but halRandGen's
    # pointers are arbitrary, so here --noDupes is only required to pick from the
    # traversed copies, and its agreement with halLiftover is reported.
    cat >"$TMP/cmp.py" <<'PYEOF'
import sys, collections
def bases(path):
    c = collections.Counter()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        s, e = int(f[1]), int(f[2])
        key = (f[0], f[5] if len(f) > 5 else ".", tuple(f[3:5]))
        if len(f) >= 12:
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            starts = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = [(s + b, s + b + l) for b, l in zip(starts, sizes)]
        else:
            blocks = [(s, e)]
        for bs, be in blocks:
            for p in range(bs, be): c[key + (p,)] += 1
    return c
mode, x, y = sys.argv[1:4]
if mode == "within":        # every base of x is a base of y
    print("ok" if not (bases(x) - bases(y)) else "OUTSIDE")
elif open(x).read() == open(y).read():
    print("bytes")
else:
    print("bases" if bases(x) == bases(y) else "DIFFERENT")
PYEOF
    H="$TMP/r7.hdf5.hal"; A="$TMP/r7.hdf5"
    genomes=($("$BIN/halStats" --genomes "$H"))
    nb=0; nbase=0; ndup=0; ndupAgree=0
    for src in "${genomes[@]:0:4}"; do
        "$BIN/halStats" --sequenceStats "$src" "$H" | awk -F', ' 'NR>1 && $2>0 {print $1"\t"$2}' >"$TMP/seqs"
        # intervals of several widths and strands, as BED6 and as BED12 with two blocks
        awk 'BEGIN{srand(7)} {n[NR]=$1; l[NR]=$2} END{for (i=0;i<40;i++){k=1+int(rand()*NR); w=1+int(rand()*2000); s=int(rand()*l[k]); if(s+w>l[k]) w=l[k]-s; if(w<1) continue;
            st=substr("+-.", 1+int(rand()*3), 1); printf "%s\t%d\t%d\tr%d\t0\t%s\n", n[k], s, s+w, i, st;
            if (w>=4) {h=int(w/2); printf "%s\t%d\t%d\tb%d\t0\t%s\t%d\t%d\t0,0,0\t2\t%d,%d,\t0,%d,\n", n[k], s, s+w, i, st, s, s+w, h-1, w-h, h >"/dev/stderr"}}}' \
            "$TMP/seqs" >"$TMP/in6.bed" 2>"$TMP/in12.bed"
        for tgt in "${genomes[@]}"; do
            [ "$src" = "$tgt" ] && continue
            for bed in in6 in12; do
                "$BIN/halLiftover" "$H" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/hal.bed" 2>/dev/null
                "$BIN/hal" liftover "$A" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/f3.bed" 2>/dev/null
                case $(python3 "$TMP/cmp.py" same "$TMP/hal.bed" "$TMP/f3.bed") in
                    bytes) nb=$((nb + 1)) ;;
                    bases) nbase=$((nbase + 1)) ;;
                    *) echo "  FAIL $src -> $tgt $bed: different aligned bases"; fails=$((fails + 1)) ;;
                esac
                "$BIN/halLiftover" "$H" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/haln.bed" --noDupes 2>/dev/null
                "$BIN/hal" liftover "$A" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/f3n.bed" --noDupes 2>/dev/null
                ndup=$((ndup + 1))
                [ "$(python3 "$TMP/cmp.py" same "$TMP/haln.bed" "$TMP/f3n.bed")" != DIFFERENT ] && ndupAgree=$((ndupAgree + 1))
                if [ "$(python3 "$TMP/cmp.py" within "$TMP/f3n.bed" "$TMP/f3.bed")" != ok ]; then
                    echo "  FAIL $src -> $tgt $bed --noDupes: maps outside the traversed copies"; fails=$((fails + 1))
                fi
            done
        done
    done
    echo "  ok   $((nb + nbase)) runs: $nb byte-identical, $nbase the same aligned bases grouped differently"
    echo "  ok   $ndup --noDupes runs pick from the traversed copies ($ndupAgree agree with halLiftover's pointers)"
fi

echo "hal maf against hal2maf:"
if [ ! -x "$BIN/hal2maf" ] || [ ! -x "$BIN/halStats" ]; then
    echo "  skip (needs hal2maf and halStats: build the maf and stats modules)"
else
    # Blocks are cut differently by design; what must agree is every column holding a
    # reference base: the reference base and every other base in it, with strand and
    # letter.  (--noDupes is left out: halRandGen's canonical-paralog pointers are
    # arbitrary, see the liftover tests.)
    cat >"$TMP/mafcols.py" <<'PYEOF'
import sys
rows = []
out = sys.stdout
def fwd(r, p):
    return p if r[2] == "+" else r[3] - 1 - p
def flush():
    if not rows: return
    anc = rows[0]
    apos, p = [], anc[1]
    for c in anc[4]:
        if c == "-": apos.append(None)
        else: apos.append(fwd(anc, p)); p += 1
    for x in apos:
        if x is not None: out.write("REF\t%s\t%d\n" % (anc[0], x))
    for r in rows[1:]:
        p = r[1]
        for i, c in enumerate(r[4]):
            if c == "-": continue
            if apos[i] is not None:
                out.write("%s\t%d\t%s\t%d\t%s\t%s\n" % (anc[0], apos[i], r[0], fwd(r, p), "+" if r[2] == anc[2] else "-", c))
            p += 1
for line in open(sys.argv[1]):
    if line.startswith("a"):
        flush(); rows = []
    elif line.startswith("s"):
        t = line.split()
        rows.append((t[1], int(t[2]), t[4], int(t[5]), t[6]))
flush()
PYEOF
    nmaf=0
    for seed in 1 7; do
        H="$TMP/r$seed.hdf5.hal"; A="$TMP/r$seed.hdf5"
        genomes=($("$BIN/halStats" --genomes "$H"))
        for ref in "${genomes[0]}" "${genomes[1]}" "${genomes[${#genomes[@]}-1]}"; do
            for opt in "" --onlyOrthologs --noAncestors; do
                if [ "$opt" = --noAncestors ] && [ "$("$BIN/halStats" "$H" | awk -F', ' -v g="$ref" '$1==g {print $2}')" != 0 ]; then
                    # an ancestral reference: both must refuse
                    if "$BIN/hal2maf" "$H" "$TMP/h.maf" --refGenome "$ref" $opt 2>/dev/null ||
                       "$BIN/hal" maf "$A" "$TMP/f.maf" --refGenome "$ref" --uncompressed $opt 2>/dev/null; then
                        echo "  FAIL seed $seed: --noAncestors with ancestral reference $ref was accepted"; fails=$((fails + 1))
                    fi
                    continue
                fi
                "$BIN/hal2maf" "$H" "$TMP/h.maf" --refGenome "$ref" $opt 2>/dev/null
                "$BIN/hal" maf "$A" "$TMP/f.maf" --refGenome "$ref" --uncompressed $opt 2>/dev/null
                python3 "$TMP/mafcols.py" "$TMP/h.maf" | LC_ALL=C sort >"$TMP/h.cols"
                python3 "$TMP/mafcols.py" "$TMP/f.maf" | LC_ALL=C sort >"$TMP/f.cols"
                if cmp -s "$TMP/h.cols" "$TMP/f.cols" && [ -s "$TMP/h.cols" ]; then
                    nmaf=$((nmaf + 1))
                else
                    echo "  FAIL seed $seed ref $ref $opt: columns differ from hal2maf"; fails=$((fails + 1))
                fi
            done
        done
    done
    echo "  ok   $nmaf exports: every reference column identical to hal2maf's, bases included"
    "$BIN/hal" maf "$TMP/r7.hdf5" "$TMP/z.maf.gz" 2>/dev/null
    "$BIN/hal" maf "$TMP/r7.hdf5" "$TMP/z.maf" --uncompressed 2>/dev/null
    check "BGZF output decompresses to the plain output" cmp "$TMP/z.maf" <(gzip -dc "$TMP/z.maf.gz")
    "$BIN/hal" maf "$TMP/r7.hdf5" "$TMP/c.maf" --uncompressed --chunk 777 2>/dev/null
    python3 "$TMP/mafcols.py" "$TMP/z.maf" | LC_ALL=C sort >"$TMP/z.cols"
    python3 "$TMP/mafcols.py" "$TMP/c.maf" | LC_ALL=C sort >"$TMP/c.cols"
    check "a small chunk size changes blocks only, not columns" cmp "$TMP/z.cols" "$TMP/c.cols"
    # a reference two levels below the root, whose chunks can be wholly parentless: the
    # walk must still reach every genome below (the total of aligned bases is the same
    # for any chunk size)
    for deep in $("$BIN/hal" stats "$TMP/r7.hdf5" --genomes); do
        p=$("$BIN/hal" stats "$TMP/r7.hdf5" --parent "$deep")
        [ -n "$p" ] && [ -n "$("$BIN/hal" stats "$TMP/r7.hdf5" --parent "$p")" ] && break
    done
    total() { "$BIN/hal" maf "$TMP/r7.hdf5" stdout --refGenome "$deep" --joinGap 0 --noDna "$@" 2>/dev/null | awk '/^s/ {t+=$4} END {print t}'; }
    check "reference $deep, chunks of 97 bases: the same aligned bases as one chunk" test "$(total)" = "$(total --chunk 97)"

    # insertion columns (--joinGap, on by default): valid MAF; no inserted base also in
    # a reference column or twice in one column; every pair of bases sharing an
    # insertion column aligned in HAL (a pair of hal2maf --global)
    cat >"$TMP/check.py" <<'PYEOF'
"""Checks for hal maf output with insertion columns.

check.py validate f.maf     MAF validity and the insertion invariants; prints counts
check.py inspairs f.maf     aligned pairs within columns that hold no reference base
check.py allpairs f.maf     aligned pairs within every column (for hal2maf --global)

A base is name:fwdpos; a pair is "a b rel" with a < b and rel + (same strand) or -.
The anchor (reference) is the first row of each block.
"""
import sys


def blocks(path):
    rows = []
    for line in open(path):
        if line.startswith("a"):
            if rows:
                yield rows
            rows = []
        elif line.startswith("s"):
            t = line.split()
            rows.append((t[1], int(t[2]), int(t[3]), t[4], int(t[5]), t[6]))
    if rows:
        yield rows


def fwd(strand, src, p):
    return p if strand == "+" else src - 1 - p


def columns(rows):
    """For each column, the list of (name, fwdpos, strand) of rows with a base there."""
    ncol = len(rows[0][5])
    cols = [[] for _ in range(ncol)]
    for name, start, size, strand, src, text in rows:
        p = start
        for i, c in enumerate(text):
            if c != "-":
                cols[i].append((name, fwd(strand, src, p), strand))
                p += 1
    return cols


def pairs(col):
    for i in range(len(col)):
        for j in range(i + 1, len(col)):
            a, b = col[i], col[j]
            if (a[0], a[1]) > (b[0], b[1]):
                a, b = b, a
            print(f"{a[0]}:{a[1]}\t{b[0]}:{b[1]}\t{'+' if a[2] == b[2] else '-'}")


mode, path = sys.argv[1], sys.argv[2]
if mode == "validate":
    bad = 0
    nblocks = inscols = insbases = refcols = 0
    inRef, inIns = set(), set()
    insDup = 0
    for rows in blocks(path):
        nblocks += 1
        ncol = len(rows[0][5])
        for name, start, size, strand, src, text in rows:
            if len(text) != ncol:
                bad += 1; print("text length differs", name, start, file=sys.stderr)
            if size != len(text) - text.count("-"):
                bad += 1; print("size != bases", name, start, size, file=sys.stderr)
            if start < 0 or start + size > src:
                bad += 1; print("outside sequence", name, start, size, src, file=sys.stderr)
        anc = rows[0][5]
        for i, col in enumerate(columns(rows)):
            if anc[i] != "-":
                refcols += 1
                for b in col[1:]:
                    inRef.add((b[0], b[1]))
            else:
                inscols += 1
                ks = [(b[0], b[1]) for b in col]
                insbases += len(ks)
                insDup += len(ks) - len(set(ks))       # twice in one column: wrong
                inIns.update(ks)
                if not col:
                    bad += 1; print("empty column", file=sys.stderr)
    both = inRef & inIns
    print(f"blocks {nblocks}, reference columns {refcols}, insertion columns {inscols}, inserted bases {insbases}")
    print(f"inserted bases also in a reference column: {len(both)}; inserted twice in one column: {insDup}; format errors: {bad}")
    sys.exit(1 if (both or insDup or bad) else 0)
elif mode in ("inspairs", "allpairs"):
    for rows in blocks(path):
        anc = rows[0][5]
        for i, col in enumerate(columns(rows)):
            if mode == "allpairs" or anc[i] == "-":
                pairs(col)
PYEOF
    # short segments, so that there are insertions short enough to join
    H="$TMP/short.hal"; A="$TMP/short"
    "$BIN/halRandGen" --testRand --seed 11 --minSegmentLength 5 --maxSegmentLength 40 --minSegments 200 \
        --maxSegments 400 "$H" >/dev/null 2>&1
    "$BIN/hal" convert "$H" "$A" --check 2000 --directory >/dev/null 2>&1
    leaf=$("$BIN/halStats" "$H" | awk -F', ' 'NR>2 && $2==0 {print $1; exit}')
    "$BIN/hal" maf "$A" "$TMP/j.maf" --refGenome "$leaf" --uncompressed 2>/dev/null
    check "joined output: valid MAF, inserted bases in no reference column, none twice in a column" \
        python3 "$TMP/check.py" validate "$TMP/j.maf"
    "$BIN/hal2maf" "$H" "$TMP/g.maf" --global 2>/dev/null
    python3 "$TMP/check.py" allpairs "$TMP/g.maf" | LC_ALL=C sort -u >"$TMP/g.pairs"
    python3 "$TMP/check.py" inspairs "$TMP/j.maf" | LC_ALL=C sort -u >"$TMP/j.pairs"
    check "insertion columns: $(wc -l <"$TMP/j.pairs") aligned pairs, all of them in hal2maf --global" \
        test "$(wc -l <"$TMP/j.pairs")" -gt 100 -a "$(LC_ALL=C comm -23 "$TMP/j.pairs" "$TMP/g.pairs" | wc -l)" -eq 0
    "$BIN/hal" maf "$A" "$TMP/j0.maf" --refGenome "$leaf" --uncompressed --joinGap 0 2>/dev/null
    python3 "$TMP/mafcols.py" "$TMP/j.maf" | LC_ALL=C sort >"$TMP/j.cols"
    python3 "$TMP/mafcols.py" "$TMP/j0.maf" | LC_ALL=C sort >"$TMP/j0.cols"
    check "joining leaves every reference column as it was" cmp "$TMP/j.cols" "$TMP/j0.cols"

    # threads, TAF, row selection -- on the short-segment alignment, which has joins
    same=0
    for ext in maf maf.gz taf taf.gz; do
        "$BIN/hal" maf "$A" "$TMP/t1.$ext" --refGenome "$leaf" --chunk 3000 2>/dev/null
        "$BIN/hal" maf "$A" "$TMP/t3.$ext" --refGenome "$leaf" --chunk 3000 --threads 3 2>/dev/null
        cmp -s "$TMP/t1.$ext" "$TMP/t3.$ext" && [ -s "$TMP/t1.$ext" ] && same=$((same + 1))
    done
    check "--threads 3: the same bytes as 1 thread, for .maf, .maf.gz, .taf and .taf.gz" test $same -eq 4
    check "the name picks the format: .maf.gz is BGZF MAF, .taf.gz BGZF TAF" \
        sh -c "gzip -dc '$TMP/t1.maf.gz' | cmp -s - '$TMP/t1.maf' && gzip -dc '$TMP/t1.taf.gz' | cmp -s - '$TMP/t1.taf' && head -1 '$TMP/t1.taf' | grep -q '^#taf'"
    # TAF read back, as taffy reads it: each ';' line starts a block whose rows are the
    # last block's, each moved past its bases, then changed by the operations in order
    cat >"$TMP/taf2maf.py" <<'PYEOF'
import sys
rle, rows, block = False, [], None
out = sys.stdout
def flush():
    if block is None: return
    out.write("a\n")
    for r, t in zip(rows, block):
        size = len(t) - t.count("-")
        out.write("s\t%s\t%d\t%d\t%s\t%d\t%s\n" % (r[0], r[1], size, r[2], r[3], t))
    out.write("\n")
for n, line in enumerate(open(sys.argv[1])):
    if n == 0:
        rle = "run_length_encode_bases:1" in line.split(); continue
    if line.startswith("#"): continue
    tok = line.split()
    k = tok.index(";") if ";" in tok else len(tok)
    bases = tok[:k]
    col = "".join(bases[i] * int(bases[i + 1]) for i in range(0, len(bases), 2)) if rle else bases[0]
    if k < len(tok):
        flush()
        if block is not None:
            for r, t in zip(rows, block): r[1] += len(t) - t.count("-")
        ops, i = tok[k + 1:], 0
        while i < len(ops):
            o, j = ops[i], int(ops[i + 1])
            if o == "d": del rows[j]; i += 2
            elif o in "is":
                r = [ops[i + 2], int(ops[i + 3]), ops[i + 4], int(ops[i + 5])]
                if o == "i": rows.insert(j, r)
                else: rows[j] = r
                i += 6
            elif o == "g": rows[j][1] += int(ops[i + 2]); i += 3
            else: sys.exit("unknown operation " + o)
        block = ["" for _ in rows]
    if len(col) != len(rows): sys.exit("line %d: %d bases for %d rows" % (n + 1, len(col), len(rows)))
    block = [t + c for t, c in zip(block, col)]
flush()
PYEOF
    "$BIN/hal" maf "$A" "$TMP/tp.taf" --refGenome "$leaf" --chunk 3000 --noRunLength 2>/dev/null
    check "TAF (run-length encoded and not) read back gives the MAF's blocks exactly" \
        sh -c "grep -v -e '^#' -e '^\$' '$TMP/t1.maf' >'$TMP/t1.body' && python3 '$TMP/taf2maf.py' '$TMP/t1.taf' | grep -v '^\$' | cmp -s - '$TMP/t1.body' && python3 '$TMP/taf2maf.py' '$TMP/tp.taf' | grep -v '^\$' | cmp -s - '$TMP/t1.body'"
    if command -v taffy >/dev/null; then
        check "taffy view -m gives the MAF back" sh -c "taffy view -m -i '$TMP/t1.taf' | grep -v -e '^#' -e '^\$' | cmp -s - '$TMP/t1.body'"
    fi
    genomes=($("$BIN/halStats" --genomes "$H"))
    other=${genomes[1]}; [ "$other" = "$leaf" ] && other=${genomes[2]}
    "$BIN/hal" maf "$A" "$TMP/x.maf" --refGenome "$leaf" --excludeGenomes "$other" --joinGap 0 2>/dev/null
    python3 "$TMP/mafcols.py" "$TMP/x.maf" | LC_ALL=C sort >"$TMP/x.cols"
    check "--excludeGenomes $other: the default columns without $other's rows" \
        cmp "$TMP/x.cols" <(awk -F'\t' -v g="$other." 'index($3, g) != 1' "$TMP/j0.cols")
    # a reference with other copies of its bases in the output
    for dup in "${genomes[@]}"; do
        "$BIN/hal" maf "$A" "$TMP/d.maf" --refGenome "$dup" --joinGap 0 2>/dev/null
        python3 "$TMP/mafcols.py" "$TMP/d.maf" | LC_ALL=C sort >"$TMP/d.cols"
        [ "$(awk -F'\t' -v g="$dup." 'index($3, g) == 1' "$TMP/d.cols" | wc -l)" -gt 0 ] && break
    done
    "$BIN/hal" maf "$A" "$TMP/s.maf" --refGenome "$dup" --refSingleCopy --joinGap 0 2>/dev/null
    python3 "$TMP/mafcols.py" "$TMP/s.maf" | LC_ALL=C sort >"$TMP/s.cols"
    check "--refSingleCopy ($dup): the default columns without $(awk -F'\t' -v g="$dup." 'index($3, g) == 1' "$TMP/d.cols" | wc -l) bases of $dup's other rows" \
        sh -c "test \$(awk -F'\t' -v g='$dup.' 'index(\$3, g) == 1' '$TMP/d.cols' | wc -l) -gt 0 && awk -F'\t' -v g='$dup.' 'index(\$3, g) != 1' '$TMP/d.cols' | cmp -s - '$TMP/s.cols'"
fi

echo "hal maf --universal:"
# the universal checker: reads a universal MAF on its own and checks every base is in
# one row at most, row 0 is + and is the backbone's (or an extra sequence), and each
# genome holds exactly the bases due: the backbone all; an ancestor its top intervals
# (all with --allAncestorRows); a leaf all but its top intervals (all with
# --leafInsertions).  Top intervals: the gaps between the runs of hal stats --runs.
#   unicheck.py check MAF GENOMES.tsv MODE      GENOMES.tsv: genome, kind, sizes file, runs file
#   unicheck.py tops SIZES RUNS                 the top intervals as BED
#   unicheck.py cols MAF...                     every column, as its sorted bases
cat >"$TMP/unicheck.py" <<'PYEOF'
import sys, collections
def intervals(sizes, runs):
    L = {}
    for l in open(sizes):
        n, v = l.split()[:2]; L[n] = int(v)
    cov = collections.defaultdict(list)
    for l in open(runs):
        f = l.split("\t")
        if len(f) >= 3: cov[f[0]].append((int(f[1]), int(f[2])))
    tops = {}
    for n, ln in L.items():
        p, g = 0, []
        for a, b in sorted(cov[n]):
            if a > p: g.append((p, a))
            p = max(p, b)
        if p < ln: g.append((p, ln))
        tops[n] = g
    return L, tops
def rows(path):
    """(block number, row number, name, fwd positions) of every row"""
    blk, k, n = -1, 0, 0
    for line in open(path):
        if line.startswith("a"): blk += 1; k = 0
        elif line.startswith("s"):
            _, name, start, size, strand, slen, text = line.split()
            start, size, slen = int(start), int(size), int(slen)
            pos, out = start, []
            for c in text:
                if c == "-": out.append(None); continue
                out.append(pos if strand == "+" else slen - 1 - pos); pos += 1
            if pos - start != size: sys.exit("%s: size %d but %d bases" % (name, size, pos - start))
            yield blk, k, name, strand, out
            k += 1
if sys.argv[1] == "tops":
    L, tops = intervals(sys.argv[2], sys.argv[3])
    for n in L:
        for a, b in tops[n]: print("%s\t%d\t%d" % (n, a, b))
elif sys.argv[1] == "cols":
    cols = collections.defaultdict(set)
    for path in sys.argv[2:]:
        for blk, k, name, strand, out in rows(path):
            for i, p in enumerate(out):
                if p is not None: cols[(path, blk, i)].add((name, p))
    for c in sorted(tuple(sorted(v)) for v in cols.values()): print(c)
else:
    maf, genomes, mode = sys.argv[2], sys.argv[3], sys.argv[4].split(",")
    kind, want, byname = {}, {}, {}
    backbone = None
    for l in open(genomes):
        g, k, sizes, runs = l.split()
        kind[g] = k
        L, tops = intervals(sizes, runs)
        for n in L: byname[g + "." + n] = (g, n, L[n])
        if k == "backbone": backbone = g
        allb = {n: [(0, L[n])] for n in L}
        if k == "backbone" or (k == "leaf" and "leaves" in mode) or (k == "anc" and "allanc" in mode): want[g] = allb
        elif k == "anc": want[g] = tops
        elif k == "excluded": want[g] = {n: [] for n in L}
        else:
            w = {}
            for n in L:
                p, g2 = 0, []
                for a, b in tops[n]:
                    if a > p: g2.append((p, a))
                    p = b
                if p < L[n]: g2.append((p, L[n]))
                w[n] = g2
            want[g] = w
    seen = collections.defaultdict(set)
    bad = []
    for blk, k, name, strand, out in rows(maf):
        g, n, ln = byname[name]
        if k == 0 and (strand != "+" or (g != backbone and kind[g] == "leaf" and "leaves" not in mode)):
            bad.append("block %d: row 0 %s %s" % (blk, name, strand))
        for p in out:
            if p is None: continue
            if (n, p) in seen[g]: bad.append("%s %d twice" % (name, p))
            seen[g].add((n, p))
    for g, w in want.items():
        exp = set((n, p) for n, iv in w.items() for a, b in iv for p in range(a, b))
        if exp != seen[g]:
            bad.append("%s: %d bases, %d due (%d missing, %d extra)" % (g, len(seen[g]), len(exp), len(exp - seen[g]), len(seen[g] - exp)))
    if bad: sys.exit("\n".join(bad[:5]))
PYEOF
"$BIN/halRandGen" --testRand --seed 11 --minSegmentLength 5 --maxSegmentLength 40 --minSegments 1000 --maxSegments 2000 \
    --maxBranchLength 3 "$TMP/u11.hal" >/dev/null 2>&1
"$BIN/hal" convert "$TMP/u11.hal" "$TMP/u11" --directory >/dev/null 2>&1
for A in "$TMP/r7.hdf5" "$TMP/u11"; do
    a=$(basename "$A")
    root=$("$BIN/hal" stats "$A" --root)
    : >"$TMP/$a.genomes"
    for g in $("$BIN/hal" stats "$A" --genomes); do
        "$BIN/hal" stats "$A" --chromSizes "$g" >"$TMP/$a.$g.sizes"
        "$BIN/hal" stats "$A" --runs "$g" >"$TMP/$a.$g.runs"
        k=anc; [ -z "$("$BIN/hal" stats "$A" --children "$g")" ] && k=leaf; [ "$g" = "$root" ] && k=backbone
        echo "$g $k $TMP/$a.$g.sizes $TMP/$a.$g.runs" >>"$TMP/$a.genomes"
    done
    # (--noDupes changes which bases are tops -- the copy rule's -- so it is checked
    # where every genome is due all its bases)
    for mode in "" --leafInsertions --allAncestorRows "--leafInsertions --allAncestorRows" "--leafInsertions --allAncestorRows --noDupes"; do
        m=$(echo "$mode" | sed 's/--leafInsertions/leaves/; s/--allAncestorRows/allanc/; s/--noDupes/nodupes/; s/ /,/g')
        "$BIN/hal" maf "$A" "$TMP/u.maf" --universal $mode --noDna --chunk 97 --threads 3 2>"$TMP/u.err"
        check "$a ${mode:-default}: every base due in one row, no other, row 0 the backbone's on +" \
            sh -c "grep -q 'every base once: checked' '$TMP/u.err' && python3 '$TMP/unicheck.py' check '$TMP/u.maf' '$TMP/$a.genomes' '$m'"
    done
    # the columns are the reference columns of the root's bases and of each ancestor's
    # top intervals, with that ancestor as the reference: hal maf's own, independently
    "$BIN/hal" maf "$A" "$TMP/u.maf" --universal --allAncestorRows --noDna 2>/dev/null
    refs=("$TMP/o.$root.maf")
    "$BIN/hal" maf "$A" "$TMP/o.$root.maf" --joinGap 0 --noDna 2>/dev/null
    while read -r g k sizes runs; do
        [ "$k" = anc ] || continue
        python3 "$TMP/unicheck.py" tops "$sizes" "$runs" >"$TMP/o.$g.bed"
        [ -s "$TMP/o.$g.bed" ] || continue
        "$BIN/hal" maf "$A" "$TMP/o.$g.maf" --refGenome "$g" --rootGenome "$g" --refTargets "$TMP/o.$g.bed" --joinGap 0 --noDna 2>/dev/null
        refs+=("$TMP/o.$g.maf")
    done <"$TMP/$a.genomes"
    check "$a: the columns are the reference columns of the root and of each ancestor's top intervals (${#refs[@]} exports)" \
        cmp <(python3 "$TMP/unicheck.py" cols "$TMP/u.maf") <(python3 "$TMP/unicheck.py" cols "${refs[@]}")
done
A="$TMP/u11"
for ext in maf maf.gz taf taf.gz; do
    "$BIN/hal" maf "$A" "$TMP/u1.$ext" --universal --leafInsertions --chunk 500 --universalBed "$TMP/u1.$ext.bed" 2>/dev/null
    "$BIN/hal" maf "$A" "$TMP/u3.$ext" --universal --leafInsertions --chunk 500 --threads 3 --universalBed "$TMP/u3.$ext.bed" 2>/dev/null
done
check "--threads 3: the same bytes as 1 thread (.maf, .maf.gz, .taf, .taf.gz, --universalBed)" \
    sh -c "for e in maf maf.gz taf taf.gz; do cmp -s '$TMP/u1.'\$e '$TMP/u3.'\$e && cmp -s '$TMP/u1.'\$e.bed '$TMP/u3.'\$e.bed && [ -s '$TMP/u1.'\$e ] || exit 1; done"
"$BIN/hal" maf "$A" "$TMP/u97.maf" --universal --leafInsertions --chunk 97 --universalBed "$TMP/u97.bed" 2>/dev/null
check "--chunk 97: the same columns, in the same order (the axis), as --chunk 500" \
    bash -c "cmp -s '$TMP/u97.bed' '$TMP/u1.maf.bed' && cmp -s <(python3 '$TMP/unicheck.py' cols '$TMP/u97.maf') <(python3 '$TMP/unicheck.py' cols '$TMP/u1.maf')"
check "headers: the universal tags and the tree" \
    sh -c "head -2 '$TMP/u1.maf' | grep -q '^##maf .*universal=flank1 leafInsertions=1 ancestorRows=own noDupes=0 chunk=500' && head -2 '$TMP/u1.maf' | grep -q '^# hal (' && head -2 '$TMP/u1.taf' | grep -q '^#taf .*universal:flank1' && head -2 '$TMP/u1.taf' | grep -q '^# hal ('"
check "TAF read back gives the MAF's blocks exactly" \
    sh -c "grep -v -e '^#' -e '^\$' '$TMP/u1.maf' >'$TMP/u1.body' && python3 '$TMP/taf2maf.py' '$TMP/u1.taf' | grep -v '^\$' | cmp -s - '$TMP/u1.body'"
if command -v taffy >/dev/null; then
    check "taffy view -m gives the MAF back" sh -c "taffy view -m -i '$TMP/u1.taf' | grep -v -e '^#' -e '^\$' | cmp -s - '$TMP/u1.body'"
fi
seq1=$("$BIN/hal" stats "$A" --sequences "$("$BIN/hal" stats "$A" --root)" | cut -d, -f2)
"$BIN/hal" maf "$A" "$TMP/us.maf" --universal --leafInsertions --chunk 500 --refSequence "$seq1" 2>/dev/null
check "--refSequence $seq1: its section of the whole export" \
    sh -c "awk -v s='$("$BIN/hal" stats "$A" --root).$seq1' '/^a/ {b = \$0 \"\\n\"; f = 1; next} /^s/ {if (f) k = \$2 == s; f = 0; b = b \$0 \"\\n\"; next} /^\$/ {if (k && b != \"\") printf \"%s\", b; b = \"\"; k = 0}' '$TMP/u1.maf' >'$TMP/us.a' && grep -v -e '^#' -e '^\$' '$TMP/us.maf' | cmp -s - '$TMP/us.a' && test -s '$TMP/us.a'"
anc=$(awk '$2 == "anc" {print $1; exit}' "$TMP/u11.genomes")
refuses "--universal with --joinGap 5" "$BIN/hal" maf "$A" "$TMP/x.maf" --universal --joinGap 5
refuses "--universal with --refGenome $anc (not the root)" "$BIN/hal" maf "$A" "$TMP/x.maf" --universal --refGenome "$anc"
refuses "--universal with --targetGenomes" "$BIN/hal" maf "$A" "$TMP/x.maf" --universal --targetGenomes "$anc"
refuses "--universal excluding an ancestor" "$BIN/hal" maf "$A" "$TMP/x.maf" --universal --excludeGenomes "$anc"

echo "hal fasta against hal2fasta:"
if [ ! -x "$BIN/hal2fasta" ] || [ ! -x "$BIN/halStats" ]; then
    echo "  skip (needs hal2fasta and halStats: build the fasta and stats modules)"
else
    H="$TMP/r7.hdf5.hal"; A="$TMP/r7.hdf5"
    genomes=($("$BIN/halStats" --genomes "$H"))
    g="${genomes[2]}"
    seq=$("$BIN/halStats" --sequenceStats "$g" "$H" | awk -F', ' 'NR==2 {print $1}')
    same() {   # same "what" options...
        local what=$1; shift
        "$BIN/hal2fasta" "$H" "$@" >"$TMP/h.fa" 2>/dev/null
        "$BIN/hal" fasta "$A" "$@" >"$TMP/f.fa" 2>/dev/null
        check "$what" cmp "$TMP/h.fa" "$TMP/f.fa"
    }
    same "a whole genome" "$g"
    same "a range of one sequence" "$g" --sequence "$seq" --start 17 --length 1234
    same "a genome range across sequences, --lineWidth 17" "$g" --start 100 --length 250000 --lineWidth 17
    same "a subtree, --ucscSequenceNames --upper" "${genomes[0]}" --subtree --ucscSequenceNames --upper
fi

echo "hal stats against halStats:"
if [ ! -x "$BIN/halStats" ]; then
    echo "  skip (needs halStats: build the stats module)"
else
    H="$TMP/r7.hdf5.hal"; A="$TMP/r7.hdf5"
    genomes=($("$BIN/halStats" --genomes "$H"))
    g="${genomes[${#genomes[@]}-1]}"; a="${genomes[1]}"
    nsame=0
    for opt in "--genomes" "--tree" "--root" "--branches" "--assembly" "--metaData" "--chromSizes $g" \
               "--sequences $a" "--bedSequences $a" "--children $a" "--parent $g" "--branchLength $g" \
               "--span $g,${genomes[2]}" "--spanRoot $g,${genomes[2]}" "--baseComp $g,1000" "--genomeMetaData $g"; do
        if [ "$("$BIN/halStats" "$H" $opt 2>&1)" = "$("$BIN/hal" stats "$A" $opt 2>&1)" ]; then
            nsame=$((nsame + 1))
        else
            echo "  FAIL hal stats $opt differs from halStats"; fails=$((fails + 1))
        fi
    done
    echo "  ok   $nsame options give halStats' output"
    check "--sequenceStats: names and lengths as halStats" \
        cmp <("$BIN/halStats" "$H" --sequenceStats "$g" | cut -d, -f1-2) <("$BIN/hal" stats "$A" --sequenceStats "$g" | cut -d, -f1-2)
fi

echo "hal coverage against hal2maf --onlyOrthologs:"
if [ ! -x "$BIN/hal2maf" ] || [ ! -x "$BIN/halStats" ]; then
    echo "  skip (needs hal2maf and halStats)"
else
    # exact coverage is the depth of each leaf's bases in hal2maf --onlyOrthologs columns
    # (liftover with paralogs); the reference covers itself once
    cat >"$TMP/covexp.py" <<'PYEOF'
import sys, collections
ref, leaves = sys.argv[1], sys.argv[2].split(",")
depth, refpos = collections.Counter(), set()
for l in sys.stdin:
    f = l.rstrip("\n").split("\t")
    if f[0] == "REF":
        refpos.add((f[1], f[2])); continue
    g = f[2].split(".")[0]
    if g in leaves: depth[(f[0], f[1], g)] += 1
hist = {g: collections.Counter() for g in leaves}
for (s, p, g), d in depth.items():
    for k in range(1, d + 1): hist[g][k] += 1
if ref in leaves: hist[ref] = collections.Counter({1: len(refpos)})
m = max([max(h) for h in hist.values() if h] or [0])
print("Genome" + "".join(", sitesCovered%dTimes" % k for k in range(1, m + 1)))
for g in leaves: print(g + "".join(", %d" % hist[g][k] for k in range(1, m + 1)))
PYEOF
    H="$TMP/r7.hdf5.hal"; A="$TMP/r7.hdf5"
    leaves=$("$BIN/halStats" "$H" | awk -F', ' 'NR>2 && NF>2 && $2==0 {print $1}')
    order=$("$BIN/hal" coverage "$A" "$("$BIN/halStats" "$H" --root)" | awk -F', ' 'NR>1 {print $1}' | paste -sd,)
    for ref in $(echo "$leaves" | head -1) $("$BIN/halStats" "$H" --root); do
        "$BIN/hal2maf" "$H" "$TMP/oo.maf" --refGenome "$ref" --onlyOrthologs 2>/dev/null
        python3 "$TMP/mafcols.py" "$TMP/oo.maf" | python3 "$TMP/covexp.py" "$ref" "$order" >"$TMP/cov.exp"
        check "exact coverage of $ref: hal2maf --onlyOrthologs column depths" cmp "$TMP/cov.exp" <("$BIN/hal" coverage "$A" "$ref")
    done
    ref=$(echo "$leaves" | head -1)
    check "--numSamples 5000: no leaf covers more than 5000 samples" \
        test "$("$BIN/hal" coverage "$A" "$ref" --numSamples 5000 --seed 3 | awk -F', ' 'NR>1 && $2>5000' | wc -l)" -eq 0
    check "--threads 3: the same tables, exact and sampled, by sequence" \
        sh -c "for o in '' '--numSamples 5000 --seed 3' '--bySequence'; do '$BIN/hal' coverage '$A' '$ref' \$o >'$TMP/c1.txt' && '$BIN/hal' coverage '$A' '$ref' \$o --threads 3 --chunk 3000 | cmp -s - '$TMP/c1.txt' || exit 1; done"
fi

echo "single file (pack):"
H="$TMP/r7.hdf5.hal"; D="$TMP/r7.hdf5"; P="$TMP/r7.pack.hal"
check "hal convert writes a single file by default, which validates" \
    sh -c "'$BIN/hal' convert '$H' '$P' --threads 3 --check 1000 >/dev/null 2>&1 && test -f '$P' && ! test -e '$P.convert-tmp' && '$BIN/hal' validate '$P'"
check "hal convert --checkOnly checks a single file" "$BIN/hal" convert "$H" "$P" --checkOnly --check 1000
"$BIN/hal" unpack "$P" "$TMP/r7.unpacked" >/dev/null 2>&1
check "unpacking the single-file conversion gives the directory conversion" diff -r "$D" "$TMP/r7.unpacked"
"$BIN/hal" pack "$D" "$TMP/r7.packed.hal" >/dev/null 2>&1
check "packing the directory gives the same file as converting to one" cmp "$P" "$TMP/r7.packed.hal"
if [ -x "$BIN/halStats" ]; then
    genomes=($("$BIN/halStats" --genomes "$H"))
    leaf=$("$BIN/halStats" "$H" | awk -F', ' 'NR>2 && NF>2 && $2==0 {print $1; exit}')
    seq=$("$BIN/halStats" --sequenceStats "$leaf" "$H" | awk -F', ' 'NR==2 {print $1}')
    printf "%s\t100\t5000\tx\t0\t+\n" "$seq" >"$TMP/p.bed"
    nsame=0
    for cmd in "stats ALN" "stats ALN --sequenceStats $leaf" "maf ALN stdout --refGenome $leaf" \
               "liftover ALN $leaf $TMP/p.bed ${genomes[1]} stdout" "fasta ALN $leaf" "coverage ALN $leaf" "validate ALN --noCrc"; do
        a=$("$BIN/hal" ${cmd//ALN/$D} 2>/dev/null | sed "s|$D|ALN|g"); b=$("$BIN/hal" ${cmd//ALN/$P} 2>/dev/null | sed "s|$P|ALN|g")
        if [ "$a" = "$b" ] && [ -n "$a" ]; then nsame=$((nsame + 1)); else echo "  FAIL hal ${cmd%% *}: single file and directory differ"; fails=$((fails + 1)); fi
    done
    echo "  ok   $nsame commands: the same output from the single file as from the directory"
fi
if [ -x "$BIN/hal3PackTest" ]; then
    cp "$P" "$TMP/crash.hal"
    check "commit protocol: abandoned update, damaged slot, torn index, interrupted create, two writers" \
        "$BIN/hal3PackTest" "$TMP/crash.hal"
fi

echo "hal update:"
P="$TMP/r7.pack.hal"
root=$("$BIN/hal" stats "$P" --root)
kids=($("$BIN/hal" stats "$P" --children "$root"))
mid="${kids[0]}"
for c in "${kids[@]}"; do [ -n "$("$BIN/hal" stats "$P" --children "$c")" ] && { mid=$c; break; }; done
grand=($("$BIN/hal" stats "$P" --children "$mid"))
prunes=""; for c in "${grand[@]}"; do prunes="$prunes --prune $c"; done
"$BIN/hal" update "$P" "$TMP/u.top.hal" $prunes >/dev/null 2>&1
"$BIN/hal" update "$P" "$TMP/u.sub.hal" --extract "$mid" >/dev/null 2>&1
"$BIN/hal" update "$TMP/u.top.hal" "$TMP/u.joined.hal" --graft "$TMP/u.sub.hal" >/dev/null 2>&1
rm -rf "$TMP/u1" "$TMP/u2"; "$BIN/hal" unpack "$P" "$TMP/u1" >/dev/null 2>&1; "$BIN/hal" unpack "$TMP/u.joined.hal" "$TMP/u2" >/dev/null 2>&1
check "split at $mid (prune below it / extract it) and graft back: every file as it was" diff -r "$TMP/u1" "$TMP/u2"
check "the extracted subtree validates" "$BIN/hal" validate "$TMP/u.sub.hal"
leaves=($("$BIN/hal" stats "$P" --genomes))
l1=""; l2=""
for g in "${leaves[@]}"; do
    if [ -z "$("$BIN/hal" stats "$P" --children "$g")" ]; then [ -z "$l1" ] && l1=$g || l2=$g; fi
done
check "--keep $l1,$l2: the spanning tree, valid" \
    sh -c "'$BIN/hal' update '$P' '$TMP/u.keep.hal' --keep $l1,$l2 >/dev/null 2>&1 && '$BIN/hal' validate '$TMP/u.keep.hal'"
cp "$P" "$TMP/u.ip.hal"
"$BIN/hal" update --inPlace "$TMP/u.ip.hal" --remove "$l2" --rename "$mid=RenamedAnc" >/dev/null 2>&1
"$BIN/hal" update "$P" "$TMP/u.new.hal" --remove "$l2" --rename "$mid=RenamedAnc" >/dev/null 2>&1
rm -rf "$TMP/u1" "$TMP/u2"; "$BIN/hal" unpack "$TMP/u.ip.hal" "$TMP/u1" >/dev/null 2>&1; "$BIN/hal" unpack "$TMP/u.new.hal" "$TMP/u2" >/dev/null 2>&1
check "--inPlace (remove + rename, one transaction) gives what a new file gets, and validates" \
    sh -c "diff -r '$TMP/u1' '$TMP/u2' && '$BIN/hal' validate '$TMP/u.ip.hal'"
check "the input is untouched by an update to a new file" cmp "$P" "$TMP/r7.packed.hal"
# a cycle of renames in place: each member must end up under its new name, not a name
# another genome has taken already
cp "$P" "$TMP/u.swap.hal"; rm -rf "$TMP/u1" "$TMP/u2"; cp -r "$D" "$TMP/u2"
swap="--rename $l1=T --rename $l2=$l1 --rename T=$l2"
"$BIN/hal" update --inPlace "$TMP/u.swap.hal" $swap >/dev/null 2>&1
"$BIN/hal" update --inPlace "$TMP/u2" $swap >/dev/null 2>&1
"$BIN/hal" unpack "$TMP/u.swap.hal" "$TMP/u1" >/dev/null 2>&1
check "--inPlace swap of two names: valid, and the single file equals the directory" \
    sh -c "'$BIN/hal' validate '$TMP/u.swap.hal' && diff -r '$TMP/u1' '$TMP/u2'"
refuses "removing an ancestor but not its children is refused" "$BIN/hal" update "$P" "$TMP/u.bad.hal" --remove "$mid"
check "... and writes nothing" test ! -e "$TMP/u.bad.hal"
refuses "grafting genomes that are already there is refused" "$BIN/hal" update "$P" "$TMP/u.bad.hal" --graft "$TMP/u.sub.hal"
refuses "replacing a genome with another alignment's root is refused" "$BIN/hal" update "$P" "$TMP/u.bad.hal" --replace "$TMP/u.sub.hal:$mid"

echo "ancestral likelihoods:"
# Independent check: in a hal maf column where every genome has at most one row, the
# column's tree is the genome tree, so Felsenstein's pruning over the leaf rows gives
# the reference's exact subtree likelihoods, and an up-and-down pass every ancestor
# row's exact posterior.  hal stores values rounded to 0.1 decade and builds each level
# from its children's stored values, so it must agree closely (KL), and its most
# probable base must be the exact one wherever that leads by more than 1.35x.
cat >"$TMP/felsenstein.py" <<'PYEOF'
import math, re, sys
mode, newick, maf = sys.argv[1], sys.argv[2], sys.argv[3]
tables = dict(zip(sys.argv[4::2], sys.argv[5::2]))     # genome -> hal fasta --likelihoods table
def parse(s):                                           # Newick -> children, branch lengths, root
    kids, parent, blen, stack, cur, pending = {}, {}, {}, [], None, []
    for t in re.findall(r"[(),;]|:[^,();]+|[^:,();]+", s):
        if t == "(": stack.append([])
        elif t == ")": pending = stack.pop()
        elif t in ",;": pass
        elif t.startswith(":"): blen[cur] = float(t[1:])
        else:
            cur = t; kids[cur] = pending
            for c in pending: parent[c] = cur
            pending = []
            if stack: stack[-1].append(cur)
    return kids, blen, [g for g in kids if g not in parent][0]
kids, blen, root = parse(open(newick).read().strip())
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
def norm(v):
    s = sum(v); return [x / s for x in v] if s > 0 else [0.25] * 4
def load(path):
    t = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.split("\t"); t[(f[0], int(f[1]))] = [float(x) for x in f[3:7]]
    return t
ours = {g: load(p) for g, p in tables.items()}
def comp(v): return [v[3], v[2], v[1], v[0]]
def msg(c, col, L):          # product over c's copies of sum_b P(b | a) L_copy(b), over a
    e = math.exp(-4 * blen[c] / 3); same, other = 0.25 + 0.75 * e, 0.25 - 0.25 * e
    if kids[c]:
        l = L[c]; s = sum(l); return [other * s + (same - other) * l[a] for a in range(4)]
    v = [1.0] * 4
    for row in col.get(c, []):                          # every copy of a leaf: evidence
        if row[0] in IDX:
            x = IDX[row[0]]; v = [v[a] * (same if a == x else other) for a in range(4)]
    return v
def up(g, col, L):
    v = [1.0] * 4
    for c in kids[g]:
        if kids[c]: up(c, col, L)
        m = msg(c, col, L); v = [v[a] * m[a] for a in range(4)]
    L[g] = norm(v)
def down(g, U, L, out, col):
    out[g] = norm([L[g][a] * U[a] for a in range(4)])
    for c in kids[g]:
        if not kids[c]: continue
        rest = list(U)
        for s in kids[g]:
            if s != c: m = msg(s, col, L); rest = [rest[a] * m[a] for a in range(4)]
        e = math.exp(-4 * blen[c] / 3); same, other = 0.25 + 0.75 * e, 0.25 - 0.25 * e
        down(c, norm([sum(rest[a] * (same if a == b else other) for a in range(4)) for b in range(4)]), L, out, col)
kl, n, flips, clear = 0.0, 0, 0, 0
def compare(exact, mine):
    global kl, n, flips, clear
    mine = norm(mine)
    kl += sum(e * (math.log2(max(e, 1e-300)) - math.log2(max(m, 1e-12))) for e, m in zip(exact, mine)); n += 1
    se = sorted(exact)
    if max(range(4), key=lambda a: exact[a]) != max(range(4), key=lambda a: mine[a]) and se[3] > se[2]:
        flips += 1
        if se[3] > 1.35 * se[2]: clear += 1
def columns(block):          # per column, genome -> its rows' (base, seq, strand, forward pos)
    for i in range(len(block[0][5])):
        col, bad = {}, False
        for name, start, size, strand, srclen, text, g, pos in block:
            if text[i] == "-": continue
            col.setdefault(g, []).append((text[i].upper(), name.split(".", 1)[1], strand, pos[i]))
            if kids.get(g) and len(col[g]) > 1: bad = True; break   # an ancestor twice: no tree
        if not bad: yield col
def flush(block):
    if not block: return
    for r in block:
        name, start, size, strand, srclen, text = r[:6]; pos, p = [], start
        for ch in text:
            pos.append(p if strand == "+" else srclen - 1 - p)
            if ch != "-": p += 1
        r.extend([name.split(".", 1)[0], pos])
    for col in columns(block):
        L = {}
        if mode == "subtree":
            g = block[0][6]; b, seq, strand, p = col[g][0]
            up(g, col, L)
            compare(L[g], ours[g][(seq, p)])
        elif root in col:
            up(root, col, L); out = {}
            down(root, [0.25] * 4, L, out, col)
            for g, rows in col.items():
                if g not in ours: continue
                b, seq, strand, p = rows[0]; mine = ours[g][(seq, p)]
                compare(out[g], mine if strand == "+" else comp(mine))
block = []
for line in open(maf):
    if line.startswith("a"): flush(block); block = []
    elif line.startswith("s"):
        f = line.split(); block.append([f[1], int(f[2]), int(f[3]), f[4], int(f[5]), f[6]])
flush(block)
print(f"{n} vectors: mean KL {kl / max(n, 1):.2e} bits, top base differs at {flips}, with a lead over 1.35x at {clear}")
sys.exit(0 if n > 200 and kl / n < 1e-2 and clear == 0 else 1)
PYEOF
# halRandGen's branches all have length 0, which would make every base certain: give
# them lengths from 0.01 to 0.3 first
bl=""; i=0
for g in $("$BIN/hal" stats "$TMP/r7.pack.hal" --genomes); do
    [ "$g" = "$("$BIN/hal" stats "$TMP/r7.pack.hal" --root)" ] && continue
    bl="$bl --branchLength $g=$(awk -v i=$i 'BEGIN {printf "%.3f", 0.01 + (i * 37 % 30) / 100}')"; i=$((i + 1))
done
"$BIN/hal" update "$TMP/r7.pack.hal" "$TMP/lb.hal" $bl >/dev/null 2>&1
P="$TMP/lb.hal"
"$BIN/hal" likelihoods "$P" "$TMP/l1.hal" >/dev/null 2>&1
"$BIN/hal" likelihoods "$P" "$TMP/l3.hal" --threads 3 >/dev/null 2>&1
check "hal likelihoods: a valid alignment, the same for any number of threads" \
    sh -c "'$BIN/hal' validate '$TMP/l1.hal' && cmp '$TMP/l1.hal' '$TMP/l3.hal'"
"$BIN/hal" stats "$P" --tree >"$TMP/tree.nwk"
ancs=$(for g in $("$BIN/hal" stats "$P" --genomes); do [ -n "$("$BIN/hal" stats "$P" --children "$g")" ] && echo "$g"; done)
tabs=""; for g in $ancs; do "$BIN/hal" fasta "$TMP/l1.hal" "$g" --likelihoods >"$TMP/$g.lik"; tabs="$tabs $g $TMP/$g.lik"; done
nanc=0
for g in $ancs; do
    "$BIN/hal" maf "$P" "$TMP/a.maf" --refGenome "$g" --rootGenome "$g" --joinGap 0 2>/dev/null
    python3 "$TMP/felsenstein.py" subtree "$TMP/tree.nwk" "$TMP/a.maf" $tabs >"$TMP/f.out" && nanc=$((nanc + 1)) || { echo "  FAIL $g: $(cat "$TMP/f.out")"; fails=$((fails + 1)); }
done
echo "  ok   subtree likelihoods of $nanc ancestors agree with Felsenstein's pruning over hal maf columns"
"$BIN/hal" posteriors "$TMP/l1.hal" "$TMP/p1.hal" >/dev/null 2>&1
check "hal posteriors: a valid alignment" "$BIN/hal" validate "$TMP/p1.hal"
tabs=""; for g in $ancs; do "$BIN/hal" fasta "$TMP/p1.hal" "$g" --likelihoods >"$TMP/$g.post"; tabs="$tabs $g $TMP/$g.post"; done
"$BIN/hal" maf "$P" "$TMP/r.maf" --joinGap 0 2>/dev/null
check "posteriors agree with an up-and-down pass over the root's hal maf columns: $(python3 "$TMP/felsenstein.py" posterior "$TMP/tree.nwk" "$TMP/r.maf" $tabs)" \
    python3 "$TMP/felsenstein.py" posterior "$TMP/tree.nwk" "$TMP/r.maf" $tabs
check "each ancestral base is its posterior's most probable (the table's base column)" \
    sh -c "cat $(for g in $ancs; do printf "'%s' " "$TMP/$g.post"; done) | awk -F'\t' '!/^#/ {m=\$4; b=\"A\"; if (\$5>m) {m=\$5; b=\"C\"} if (\$6>m) {m=\$6; b=\"G\"} if (\$7>m) {m=\$7; b=\"T\"} if (toupper(\$3) != b && m > 0.26) {v=0; for (i=4;i<=7;i++) if (\$i==m) v++; if (v==1) bad++}} END {exit bad>0}'"
refuses "hal posteriors needs subtree likelihoods first" "$BIN/hal" posteriors "$P" "$TMP/p.bad.hal"
refuses "... and posteriors are not subtree likelihoods" "$BIN/hal" posteriors "$TMP/p1.hal" "$TMP/p.bad.hal"
# the likelihoods go with their genomes: split, graft back, every file as it was
"$BIN/hal" update "$TMP/l1.hal" "$TMP/lt.hal" $prunes >/dev/null 2>&1
"$BIN/hal" update "$TMP/l1.hal" "$TMP/ls.hal" --extract "$mid" >/dev/null 2>&1
"$BIN/hal" update "$TMP/lt.hal" "$TMP/lj.hal" --graft "$TMP/ls.hal" >/dev/null 2>&1
rm -rf "$TMP/u1" "$TMP/u2"; "$BIN/hal" unpack "$TMP/l1.hal" "$TMP/u1" >/dev/null 2>&1; "$BIN/hal" unpack "$TMP/lj.hal" "$TMP/u2" >/dev/null 2>&1
check "split and graft back with likelihoods: every file as it was" diff -r "$TMP/u1" "$TMP/u2"
check "pruning below an ancestor notes that its subtree likelihoods are out of date" \
    sh -c "'$BIN/hal' update '$TMP/l1.hal' '$TMP/ln.hal' $prunes 2>&1 | grep -q 'subtree likelihoods of $mid are out of date'"
cp "$TMP/l1.hal" "$TMP/lr.hal"
check "renaming in place renames the likelihoods too (no copy)" \
    sh -c "'$BIN/hal' update --inPlace '$TMP/lr.hal' --rename $mid=RenamedAnc >/dev/null 2>&1 && '$BIN/hal' validate '$TMP/lr.hal' >/dev/null && '$BIN/hal' ls '$TMP/lr.hal' | grep -q 'RenamedAnc.hlk'"
rm -rf "$TMP/ld"; "$BIN/hal" unpack "$TMP/l1.hal" "$TMP/ld" >/dev/null 2>&1
lik=$(ls "$TMP/ld/genomes/"*.hlk | head -1)
printf '\x55' | dd of="$lik" bs=1 seek=200 conv=notrunc 2>/dev/null
refuses "hal validate catches a damaged likelihood file (even with --noCrc)" "$BIN/hal" validate "$TMP/ld" --noCrc
rm "$lik"
refuses "... and a missing one" "$BIN/hal" validate "$TMP/ld"

echo "the public API (inc/hal3/):"
# a program that includes only hal3/hal3.h, built with no liftmap, twobit64 or HAL
# header on the include path (so the public headers expose none) and linked without
# libHal or HDF5
ROOT=$(cd "$BIN/.." && pwd)
cat >"$TMP/api.cpp" <<'CPPEOF'
#include "hal3/hal3.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
using namespace hal3::api;
int run(char **argv) {
    auto aln = Alignment::open(argv[1]);
    std::string what = argv[2];
    if (what == "genomes") {   /* name, parent, length, sequences; every sequence's table entry */
        for (const Genome *g : aln->genomes()) {
            printf("%s\t%s\t%lld\t%lld\n", g->name().c_str(), g->parent() ? g->parent()->name().c_str() : "-",
                   (long long)g->length(), (long long)g->numSequences());
            int64_t at = 0;
            for (uint32_t i = 0; i < (uint32_t)g->numSequences(); i++) {
                if (g->sequenceStart(i) != at || g->findSequence(g->sequenceName(i)) != (int64_t)i) return 1;
                if (g->sequenceLength(i) > 0 && (g->sequenceAt(at) != i || g->sequenceAt(at + g->sequenceLength(i) - 1) != i)) return 1;
                at += g->sequenceLength(i);
            }
            if (at != g->length()) return 1;
        }
        return aln->root() == aln->genomes()[0] && aln->genome("no such genome") == nullptr ? 0 : 1;
    }
    const Genome *src = aln->genome(argv[3]), *tgt = aln->genome(argv[4]);
    if (what == "lift") {      /* BED on stdin, one interval per line, lifted: BED3 of the pieces */
        Lifter lifter(aln);
        std::vector<LiftPiece> pieces;
        char name[1024];
        long long s, e;
        while (scanf("%1023s %lld %lld%*[^\n]", name, &s, &e) == 3) {
            lifter.lift(src, (uint32_t)src->findSequence(name), s, e, tgt, strcmp(argv[5], "1") == 0, pieces);
            for (const LiftPiece &p : pieces)
                printf("%s\t%lld\t%lld\n", tgt->sequenceName(p.tgtSeq).c_str(), (long long)p.tgtStart, (long long)(p.tgtStart + p.length));
        }
        return 0;
    }
    if (what == "walk") {      /* bases of every printed genome for [s, e) of sequence argv[4] of src */
        ReferenceWalk walk(aln, src);
        std::map<const Genome *, std::vector<LiftPiece>> pieces;
        walk.lift((uint32_t)src->findSequence(argv[4]), {{atoll(argv[5]), atoll(argv[6])}}, pieces);
        long long n = 0;
        for (auto &gp : pieces)
            for (const LiftPiece &p : gp.second) n += p.length;
        printf("%lld\n", n);
        return 0;
    }
    return 2;
}
int main(int, char **argv) {
    try {
        return run(argv);
    } catch (Error &e) {
        fprintf(stderr, "error (kind %d): %s\n", (int)e.kind(), e.what());
    } catch (std::exception &e) {
        fprintf(stderr, "error: %s\n", e.what());
    }
    return 1;
}
CPPEOF
curl=$(nm "$ROOT/../liftmap/libliftmap.a" 2>/dev/null | grep -q ' U curl_easy_init' && echo -lcurl)
check "a program using only hal3/hal3.h builds, without liftmap, twobit64 or HAL headers or libHal" \
    ${CXX:-g++} -std=c++11 -O1 -I"$ROOT/format3/inc" "$TMP/api.cpp" -o "$TMP/api" \
    "$ROOT/lib/libHal3.a" "$ROOT/../twobit64/libtwobit64.a" "$ROOT/../liftmap/libliftmap.a" -lz -lpthread $curl
A="$TMP/r7.hdf5"
check "genomes, tree, lengths and sequence tables as hal stats" \
    sh -c "'$TMP/api' '$A' genomes >'$TMP/api.g' && '$BIN/hal' stats '$A' | sed -n '/^GenomeName/,\$p' | tail -n +2 | awk -F', ' 'NF > 1 {print \$1\"\t\"\$3\"\t\"\$4}' >'$TMP/api.s' && cut -f1,3,4 '$TMP/api.g' | cmp -s - '$TMP/api.s'"
cat >"$TMP/bases.py" <<'PYEOF'
import sys, collections
def bases(path):
    c = collections.Counter()
    for line in open(path):
        f = line.split()
        for p in range(int(f[1]), int(f[2])): c[(f[0], p)] += 1
    return c
print("same" if bases(sys.argv[1]) == bases(sys.argv[2]) else "DIFFERENT")
PYEOF
genomes=($("$BIN/hal" stats "$A" --genomes))
nl=0; bad=0
for src in "${genomes[@]:1:3}"; do
    "$BIN/hal" stats "$A" --chromSizes "$src" | awk 'BEGIN{srand(5)} $2>100 {for (i=0;i<5;i++){s=int(rand()*($2-100)); printf "%s\t%d\t%d\n", $1, s, s+1+int(rand()*99)}}' >"$TMP/api.bed"
    for tgt in "${genomes[@]}"; do
        [ "$src" = "$tgt" ] && continue
        for d in 1 0; do
            "$TMP/api" "$A" lift "$src" "$tgt" $d <"$TMP/api.bed" >"$TMP/api.out"
            "$BIN/hal" liftover "$A" "$src" "$TMP/api.bed" "$tgt" "$TMP/api.ref" $([ $d = 0 ] && echo --noDupes) 2>/dev/null
            nl=$((nl + 1)); [ "$(python3 "$TMP/bases.py" "$TMP/api.out" "$TMP/api.ref")" = same ] || bad=$((bad + 1))
        done
    done
done
check "Lifter: the bases hal liftover gives, $nl source/target/copy-rule runs" test $bad -eq 0
leaf=$("$BIN/hal" stats "$A" | sed -n '/^GenomeName/,$p' | awk -F', ' '$2==0 {print $1; exit}')
seq=$("$BIN/hal" stats "$A" --chromSizes "$leaf" | sort -k2,2nr | head -1 | cut -f1)
check "ReferenceWalk: as many bases as hal maf --joinGap 0 aligns ($leaf $seq:100-2100)" \
    sh -c "test \"\$('$TMP/api' '$A' walk '$leaf' '$seq' 100 2100)\" = \"\$('$BIN/hal' maf '$A' /dev/null --refGenome '$leaf' --refSequence '$seq' --start 100 --length 2000 --joinGap 0 2>&1 | sed -n 's/.* rows, \([0-9]*\) aligned bases.*/\1/p')\""
refuses "Lifter refuses an interval outside the sequence" "$TMP/api" "$A" lift "$leaf" "$leaf" 1 <<<"$seq 0 999999999"

if [ $fails -eq 0 ]; then echo "ALL PASS (0 failures)"; else echo "FAILURES ($fails)"; fi
exit $((fails > 0))
