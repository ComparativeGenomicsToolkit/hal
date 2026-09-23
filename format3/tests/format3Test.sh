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

if [ $fails -eq 0 ]; then echo "ALL PASS (0 failures)"; else echo "FAILURES ($fails)"; fi
exit $((fails > 0))
