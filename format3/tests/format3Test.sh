#!/bin/bash
# HAL format 3 tests: generated alignments only (halRandGen), so they need nothing
# outside this repository and its sibling liftmap / twobit64 checkouts.
#
#   - halRandGen alignments (hdf5 and mmap back ends, several seeds) convert, pass
#     halToFormat3's own check against the input, and pass hal3Validate
#   - the same alignment gives byte-identical genome files from either back end
#   - an existing alignment is not overwritten
#   - hal3Validate catches corruption, missing and orphaned files, a swapped genome
#     file, a wrong parent, and a mis-ordered roster
#   - halToFormat3 --checkOnly catches a genome's DNA replaced by another's
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
            "$BIN/halToFormat3" "$TMP/r$seed.$fmt.hal" "$TMP/r$seed.$fmt" --threads 3 --check 2000
        check "seed $seed, $fmt: hal3Validate passes" "$BIN/hal3Validate" "$TMP/r$seed.$fmt"
    done
    check "seed $seed: hdf5 and mmap inputs give identical genome files" \
        diff -r "$TMP/r$seed.hdf5/genomes" "$TMP/r$seed.mmap/genomes"
done
refuses "an existing alignment is not overwritten" "$BIN/halToFormat3" "$TMP/r1.hdf5.hal" "$TMP/r1.hdf5"
check "--checkOnly passes on an intact alignment" "$BIN/halToFormat3" "$TMP/r1.hdf5.hal" "$TMP/r1.hdf5" --checkOnly

echo "hal3Validate catches:"
A="$TMP/r7.hdf5"
G=($(ls "$A/genomes" | sed -n 's/\.hgn$//p'))       # genome file stems
fresh() { rm -rf "$TMP/t"; cp -r "$A" "$TMP/t"; }
fresh; printf '\x55' | dd of="$TMP/t/genomes/${G[1]}.hgn" bs=1 seek=300 conv=notrunc 2>/dev/null
refuses "a corrupted .hgn" "$BIN/hal3Validate" "$TMP/t"
fresh; printf '\x55' | dd of="$TMP/t/genomes/${G[1]}.hgn" bs=1 seek=300 conv=notrunc 2>/dev/null
refuses "a corrupted .hgn, even with --noCrc" "$BIN/hal3Validate" "$TMP/t" --noCrc
fresh; rm "$TMP/t/genomes/${G[2]}.2bit"
refuses "a missing .2bit" "$BIN/hal3Validate" "$TMP/t"
fresh; cp "$TMP/t/genomes/${G[1]}.hgn" "$TMP/t/genomes/extra.hgn"
refuses "an orphaned genome file" "$BIN/hal3Validate" "$TMP/t"
fresh; cp "$TMP/t/genomes/${G[1]}.hgn" "$TMP/t/genomes/${G[2]}.hgn"
refuses "one genome's .hgn in place of another's (--noCrc)" "$BIN/hal3Validate" "$TMP/t" --noCrc
# re-parent the first non-root genome (a child of the root) under the last one listed
fresh; awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {if ($1=="genome") last=$2; next}
    $1=="genome" && $3!="-" && !done {$3=last; done=1} {print}' "$A/manifest" "$A/manifest" >"$TMP/t/manifest"
refuses "a wrong parent in the roster" "$BIN/hal3Validate" "$TMP/t"
fresh; awk 'NR==1{print; next} /^genome\t/ && !held {held=$0; next} {print} END{print held}' "$A/manifest" >"$TMP/t/manifest"
refuses "the root not listed first" "$BIN/hal3Validate" "$TMP/t"

echo "halToFormat3 --checkOnly catches:"
fresh; cp "$TMP/t/genomes/${G[1]}.2bit" "$TMP/t/genomes/${G[2]}.2bit"
refuses "one genome's DNA in place of another's" "$BIN/halToFormat3" "$TMP/r7.hdf5.hal" "$TMP/t" --checkOnly

echo "hal3Liftover against halLiftover:"
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
                "$BIN/hal3Liftover" "$A" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/f3.bed" 2>/dev/null
                case $(python3 "$TMP/cmp.py" same "$TMP/hal.bed" "$TMP/f3.bed") in
                    bytes) nb=$((nb + 1)) ;;
                    bases) nbase=$((nbase + 1)) ;;
                    *) echo "  FAIL $src -> $tgt $bed: different aligned bases"; fails=$((fails + 1)) ;;
                esac
                "$BIN/halLiftover" "$H" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/haln.bed" --noDupes 2>/dev/null
                "$BIN/hal3Liftover" "$A" "$src" "$TMP/$bed.bed" "$tgt" "$TMP/f3n.bed" --noDupes 2>/dev/null
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

if [ $fails -eq 0 ]; then echo "ALL PASS (0 failures)"; else echo "FAILURES ($fails)"; fi
exit $((fails > 0))
