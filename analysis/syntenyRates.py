#!/usr/bin/env python3
# Gene synteny rate script.
# Fails on overlapping BED entries; and possibly on entries that don't
# preserve their structure in the target. (i.e. don't preserve order
# and orientation)
import sys
import os
import argparse
import itertools
from operator import itemgetter
from sonLib import bioio

# Useful itertools recipe
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def uniq(l):
    """uniquify a list of lists (can't use set since list is
    unhashable)"""
    for i, item in enumerate(l):
        if item not in l[:i]:
            yield item

def sortBed(bedFile):
    """Sort bed entries by sequence/start. Discard overlapping entries.
    Returns bed lines."""
    entries = [x.split() for x in [x for x in bedFile.read().split("\n") if x != ""]]
    if len(entries[0]) < 6:
        raise RuntimeError("BED file must have strand information")
    # stable sort
    entries.sort(key=lambda x: int(x[1]))
    entries.sort(key=lambda x: x[0])
    # Check for overlaps
    toDelete = []
    # needed so that overlaps can be detected correctly.
    prevValidEnd = int(entries[0][2])
    for i, j in pairwise(entries):
        if i[0] != j[0]:
            prevValidEnd = int(j[1])
            continue
        if int(i[2]) >= int(j[1]) or prevValidEnd > int(j[1]):
            # can't (or rather won't) make this an assert since bed12
            # entries can overlap without their exons
            # overlapping. Still we'll discard them as they don't fit
            # into our synteny definition.
            sys.stderr.write("WARNING: discarding overlapping lines %s,\n"
                             "%s\n" % ("\t".join(i), "\t".join(j)))
            toDelete.append(i)
            toDelete.append(j)
        else:
            prevValidEnd = int(i[2])
    for i in uniq(toDelete):
        entries.remove(i)
    return entries

def liftover(bedLine, opts):
    """returns dict of bed lines keyed by sequence"""
    srcFile = bioio.getTempFile("syntenyRatesSrc.bed")
    open(srcFile, 'w').write("\t".join(map(str, bedLine)) + '\n')
    outFile = bioio.getTempFile("syntenyRatesDest.bed")
    cmd = "halLiftover %s %s %s %s %s" % (opts.halFile, opts.srcGenome,
                                          srcFile, opts.destGenome, outFile)
    bioio.system(cmd)
    outBedLines = [x.split() for x in [x for x in open(outFile).read().split("\n") if x != ""]]
    outBedDict = {}
    if len(outBedLines) > 1:
        # a duplication or folding up
        chrs = list(map(itemgetter(0), outBedLines))
        for chr in set(chrs):
            lines = [x for x in outBedLines if x[0] == chr]
            strands = list(map(itemgetter(5), lines))
            if len(set(strands)) != 1:
                # maps to two different strands, ignore this chr
                print("POSSIBLYBAD: maps to two different strands on chr")
            elif opts.mergeBedLines:
                # merge bed lines into 1 larger bed line enclosing them
                minStart = min([int(x[1]) for x in lines])
                maxEnd = max([int(x[2]) for x in lines])
                print("MERGED on chr %s -- distance %d" % (chr, maxEnd - minStart))
                outBedDict[chr] = [lines[0][0], minStart, maxEnd,
                                   lines[0][3], lines[0][4],
                                   lines[0][5]]
    elif len(outBedLines) == 0:
        # no map to target at all
        print("INVALID: no map to target")
    else:
        # 1 bed line
        outBedDict[outBedLines[0][0]] = outBedLines[0]
    os.remove(srcFile)
    os.remove(outFile)
    return outBedDict

def compareLines(i, iLifted, j, jLifted):
    """Take two pairs of lines (guaranteed to be on same chr, and sorted
    by increasing start on the positive strand of the source genome)
    and return true if they don't break synteny.
    """
    # FIXME: do this somewhere else
    i[1] = int(i[1])
    i[2] = int(i[2])
    iLifted[1] = int(iLifted[1])
    iLifted[2] = int(iLifted[2])
    j[1] = int(j[1])
    j[2] = int(j[2])
    jLifted[1] = int(jLifted[1])
    jLifted[2] = int(jLifted[2])

    assert(i[0] == j[0])
    assert(iLifted[0] == jLifted[0])
    assert(i[1] <= j[1])
    sameOrientation = i[5] == j[5]
    sameOrientationInTarget = iLifted[5] == jLifted[5]
    if sameOrientation != sameOrientationInTarget:
        # I.e. the relative orientations should be the same in source
        # and target. If not there has been a synteny break.
        print("BREAK: Relative orientations not equal in source v. target")
        return False
    # Check for any overlap
    if iLifted[2] >= jLifted[1] and iLifted[2] <= jLifted[2] or \
       iLifted[1] >= jLifted[1] and iLifted[1] <= jLifted[2] or \
       iLifted[1] <= jLifted[1] and iLifted[2] >= jLifted[2]:
        print("BREAK: Overlap detected")
        return False
    # Check for ordering
    inverted = i[5] != iLifted[5]
    ordered = iLifted[1] > jLifted[1] if inverted else iLifted[1] < jLifted[1]
    if not ordered:
        print("BREAK: Order changed")
    return ordered

def getNumSyntenies(bedLines, opts):
    """Takes a list of sorted & split bed lines, returns (num syntenies,
    num valid pairs).
    """
    numValidPairs = 0
    numSyntenies = 0
    for i, j in pairwise(bedLines):
        if i[0] != j[0]:
            print("INVALID: on different query chr")
            # different chr, ignore
            continue
        iLiftedDict = liftover(i, opts)
        jLiftedDict = liftover(j, opts)
        # get intersection of chrs lifted to
        chrs = [x for x in iLiftedDict if x in jLiftedDict]
        if len(chrs) > 1:
            # unable to resolve liftover, skip
            print("INVALID: multiple possible target chrs")
            continue
        elif len(chrs) == 0:
            # Maps to different chr, ignore.
            print("INVALID: on different target chr")
            continue
        # There's just one target chr that both map to.
        iLifted = iLiftedDict[chrs[0]]
        jLifted = jLiftedDict[chrs[0]]
        numValidPairs += 1
        if compareLines(i, iLifted, j, jLifted):
            numSyntenies += 1
        else:
            print("synteny break:")
            print(i)
            print(iLifted)
            print(j)
            print(jLifted)
    return (numSyntenies, numValidPairs)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("halFile", help="HAL alignment")
    parser.add_argument("srcGenome", help="Query genome")
    parser.add_argument("bedFile", help="Bed file of genes on query genome")
    parser.add_argument("destGenome", help="Target genome")
    parser.add_argument("--mergeBedLines", action='store_true',
                        help="Merge lifted-over bed lines into a larger bed "
                        "line containing them", default=False)
    opts = parser.parse_args()
    bedLines = sortBed(open(opts.bedFile))
    (numSyntenies, numValidPairs) = getNumSyntenies(bedLines, opts)
    print("gene pair synteny rate: %f, num syntenies: %d, num pairs: %d" % (float(numSyntenies)/numValidPairs, numSyntenies, numValidPairs))
    return 0

if __name__ == '__main__':
    sys.exit(main())
