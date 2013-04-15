#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Benchmark halLodExtract by extracting with a bunch of different
stepsizes then seeing how much smaller they get.  Also (optionally)
use mafComparator (https://github.com/dentearl/mafTools) together with
comparatorSummarizer.py (https://github.com/dentearl/mwgAlignAnalysis)
"""
import argparse
import os
import sys
import copy
import subprocess
import time
from collections import defaultdict

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments

def runHalLodExtract(inHalPath, outHalPath, step):
    runShellCommand("halLodExtract %s %s %s" % (inHalPath, outHalPath, step))

def getHalTotalSegments(halPath):
    total = (0, 0)
    for genome in getHalGenomes(halPath):
        numSegs = getHalNumSegments(halPath, genome)
        total = (total[0] + numSegs[0], total[1] + numSegs[1])
    return total

def makePath(inHalPath, outDir, step, name, ext):
    inFileName = os.path.splitext(os.path.basename(inHalPath))[0]
    outPath = os.path.join(outDir, "%s_%s_%i.%s" % (inFileName,
                                                    name, step, ext))
    return outPath

def makeMaf(inHalPath, outDir, step, overwrite, doMaf):
    srcHalPath = inHalPath
    if step > 0:
        srcHalPath = makePath(inHalPath, outDir, step, "lod", "hal")
    outMafPath = makePath(inHalPath, outDir, step, "out", "maf")
    if doMaf and (overwrite or not os.path.isfile(outMafPath)):
        runShellCommand("hal2maf %s %s --maxRefGap 100" % (srcHalPath,
                                                           outMafPath))

def compMaf(inHalPath, outDir, step, overwrite, doMaf):
    srcMaf = makePath(inHalPath, outDir, 0, "out", "maf")
    tgtMaf = makePath(inHalPath, outDir, step, "out", "maf")
    xmlPath = makePath(inHalPath, outDir, step, "comp", "xml")
    sumPath = makePath(inHalPath, outDir, step, "comp", "txt")
    if doMaf and (overwrite or not os.path.isfile(xmlPath)):
        runShellCommand("mafComparator --maf1 %s --maf2 %s --out %s --samples 100000" % (
            srcMaf, tgtMaf, xmlPath))
        runShellCommand("comparatorSummarizer.py --xml %s > %s " % (xmlPath,
                                                                    sumPath))
    xmlNearPath = makePath(inHalPath, outDir, step, "comp_near", "xml")
    sumNearPath = makePath(inHalPath, outDir, step, "comp_near", "txt")
    if doMaf and (overwrite or not os.path.isfile(xmlNearPath)):
        runShellCommand(
            "mafComparator --maf1 %s --maf2 %s --out %s --near %d --samples 100000" % (
                srcMaf, tgtMaf, xmlNearPath, int(step)))
        runShellCommand("comparatorSummarizer.py --xml %s > %s " % (
            xmlNearPath, sumNearPath))


def getPrecisionRecall(inHalPath, outDir, step, doMaf):
    if doMaf:
        sumPath = makePath(inHalPath, outDir, step, "comp", "txt")
        sumFile = open(sumPath, "r")
        line = sumFile.next()
        line = sumFile.next()
        line = sumFile.next()
        tokens = line.split()
        assert tokens[2] == "self)"

        sumNearPath = makePath(inHalPath, outDir, step, "comp_near", "txt")
        sumNearFile = open(sumNearPath, "r")
        line = sumNearFile.next()
        line = sumNearFile.next()
        line = sumNearFile.next()
        tokensNear = line.split()
        assert tokensNear[2] == "self)"

        return [float(tokens[3]), float(tokens[4]),
                float(tokensNear[3]), float(tokensNear[4])]
    elif step == 0:
        return [1., 1., 1., 1.]
    else:
        return [0., 0., 0., 0.]

def getScanTime(inHalPath, outDir, step):
    srcHalPath = inHalPath
    if step > 0:
        srcHalPath = makePath(inHalPath, outDir, step, "lod", "hal")
    genomes = getHalGenomes(inHalPath)
    assert len(genomes) > 1
    genName = genomes[1]
    bedPath = makePath(inHalPath, outDir, step, genName, "bed")
    print (step, bedPath)
    t1 = time.time()
    runShellCommand("halBranchMutations %s %s --refFile %s" % (
        srcHalPath, genName, bedPath))
    elapsedTime = time.time() - t1
    return [elapsedTime]
    
def printTable(table):
    print "Step, kb,, nTop,, nBottom,, Prec., Recall, PrecNear, RecallNear, ScanTime"
    for step, data in sorted(table.items()):
        line = "%d" % step
        idx = 0
        for elem in data:
            line += ", %s" % str(elem)
            if idx <= 2:
                orig = table[0][idx]
                frac = 1
                if orig > 0:
                    frac = float(elem) / float(orig)
                line += ", %f" % frac
            idx += 1
        print line
    
def runSteps(inHalPath, outDir, steps, overwrite, doMaf):
    table = defaultdict(list)
    makeMaf(inHalPath, outDir, 0, overwrite, doMaf)

    table[0] = [os.path.getsize(inHalPath) / 1024]
    table[0] += list(getHalTotalSegments(inHalPath))
    table[0] += getPrecisionRecall(inHalPath, outDir, 0, False)
    table[0] += getScanTime(inHalPath, outDir, 0)

    for step in steps:    
        outPath = makePath(inHalPath, outDir, step, "lod", "hal")
        
        if overwrite is True or not os.path.isfile(outPath):
            runHalLodExtract(inHalPath, outPath, step)

        makeMaf(inHalPath, outDir, step, overwrite, doMaf)
        compMaf(inHalPath, outDir, step, overwrite, doMaf)
        
        table[step] = [os.path.getsize(outPath) / 1024]
        table[step] += list(getHalTotalSegments(outPath))
        table[step] += getPrecisionRecall(inHalPath, outDir, step, doMaf)
        table[step] += getScanTime(inHalPath, outDir, step)

    return table

    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("hal", help="input hal")
    parser.add_argument("outDir", help="output dir")
    parser.add_argument("steps", help="comma-separated list of stepsizes")
    parser.add_argument("--overwrite",action="store_true", default=False)
    parser.add_argument("--maf",action="store_true", default=False)
        
    args = parser.parse_args()

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    steps = [int(x) for x in args.steps.split(",")]
    table = runSteps(args.hal, args.outDir, steps, args.overwrite, args.maf)
#    print table
    printTable(table)
if __name__ == "__main__":
    sys.exit(main())
