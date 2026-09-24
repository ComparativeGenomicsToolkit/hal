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
refuses() {   # refuses "what" command...  (must exit non-zero)
    local what=$1; shift
    if "$@" >"$TMP/out" 2>&1; then echo "  FAIL $what (was accepted)"; fails=$((fails + 1)); else echo "  ok   $what"; fi
}

echo "conversion:"
for seed in 1 7 23; do
    for fmt in hdf5 mmap; do
        "$BIN/halRandGen" --testRand --seed $seed --format $fmt "$TMP/r$seed.$fmt.hal" >/dev/null 2>&1
        check "seed $seed, $fmt: converts and passes the check against HAL" \
            "$BIN/hal" convert "$TMP/r$seed.$fmt.hal" "$TMP/r$seed.$fmt" --threads 3 --check 2000
        check "seed $seed, $fmt: hal3Validate passes" "$BIN/hal" validate "$TMP/r$seed.$fmt"
    done
    check "seed $seed: hdf5 and mmap inputs give identical genome files" \
        diff -r "$TMP/r$seed.hdf5/genomes" "$TMP/r$seed.mmap/genomes"
done
refuses "an existing alignment is not overwritten" "$BIN/hal" convert "$TMP/r1.hdf5.hal" "$TMP/r1.hdf5"
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
    "$BIN/hal" convert "$H" "$A" --check 2000 >/dev/null 2>&1
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
fi

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
fi

if [ $fails -eq 0 ]; then echo "ALL PASS (0 failures)"; else echo "FAILURES ($fails)"; fi
exit $((fails > 0))
