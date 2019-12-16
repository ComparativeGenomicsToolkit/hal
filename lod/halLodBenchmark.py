#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

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

from hal.lod.halLodInterpolate import runHalLodExtract
from hal.lod.halLodInterpolate import makePath
from hal.lod.halLodInterpolate import getSteps

def getHalTotalSegments(halPath):
    total = (0, 0)
    for genome in getHalGenomes(halPath):
        numSegs = getHalNumSegments(halPath, genome)
        total = (total[0] + numSegs[0], total[1] + numSegs[1])
    return total

def makeMaf(inHalPath, outDir, step, overwrite, doMaf):
    srcHalPath = inHalPath
    if step > 0:
        srcHalPath = makePath(inHalPath, outDir, step, "lod", "hal")
    outMafPath = makePath(inHalPath, outDir, step, "out", "maf")
    if doMaf and (overwrite or not os.path.isfile(outMafPath)):
        runShellCommand("hal2maf %s %s" % (srcHalPath, outMafPath))

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
        line = next(sumFile)
        line = next(sumFile)
        line = next(sumFile)
        tokens = line.split()
        assert tokens[2] == "self)"

        sumNearPath = makePath(inHalPath, outDir, step, "comp_near", "txt")
        sumNearFile = open(sumNearPath, "r")
        line = next(sumNearFile)
        line = next(sumNearFile)
        line = next(sumNearFile)
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
    t1 = time.time()
    runShellCommand("halBranchMutations %s %s --refFile %s" % (
        srcHalPath, genName, bedPath))
    elapsedTime = time.time() - t1
    return [elapsedTime]
    
def printTable(table):
    print("Step, kb,, nTop,, nBottom,, Prec., Recall, PrecNear, RecallNear, ScanTime")
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
        print(line)
    
def runSteps(inHalPath, outDir, maxBlock, scale, steps, overwrite, doMaf,
             keepSeq, trans, inMemory):
    table = defaultdict(list)
    makeMaf(inHalPath, outDir, 0, overwrite, doMaf)

    table[0] = [os.path.getsize(inHalPath) / 1024]
    table[0] += list(getHalTotalSegments(inHalPath))
    table[0] += getPrecisionRecall(inHalPath, outDir, 0, False)
    table[0] += getScanTime(inHalPath, outDir, 0)

    if steps is None:
        steps =  getSteps(inHalPath, maxBlock, scale)
    for stepIdx in range(1,len(steps)):
        step = steps[stepIdx]
        outPath = makePath(inHalPath, outDir, step, "lod", "hal")
        
        srcPath = inHalPath
        if trans is True and stepIdx > 1:
            srcPath = makePath(inHalPath, outDir,  steps[stepIdx-1],
                               "lod", "hal")        
        
        if overwrite is True or not os.path.isfile(outPath):
            stepScale = (scale ** stepIdx)
            runHalLodExtract(srcPath, outPath, stepScale, keepSeq, inMemory)

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
    parser.add_argument("--maxBlock",
                        help="maximum desired number of blocks to ever " 
                        "display at once.", type=int,
                        default=500)
    parser.add_argument("--scale",
                        help="scaling factor between two successive levels"
                        " of detail", type=float,
                        default=10.0)
    parser.add_argument("--steps",
                        help="comma-separated list of sampling steps to test"
                        " for each level of "
                        "detail.  Overrides --scale and --maxBlock options",
                        type=str, default=None)
    parser.add_argument("--overwrite",action="store_true", default=False)
    parser.add_argument("--maf",action="store_true", default=False)
    parser.add_argument("--keepSequences",action="store_true", default=False)
    parser.add_argument("--trans", help="Generate level of detail X from "
                        "X-1.  By default, all levels of detail are generated "
                        "from the original HAL (X=0)",
                        action="store_true", default=False)
    parser.add_argument("--inMemory", help="Load entire hdf5 arrays into "
                        "memory, overriding cache.",
                        action="store_true", default=False)

        
    args = parser.parse_args()

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    if args.maf is True:
        args.keepSequences = True

    steps = None
    if args.steps is not None:
        steps = [0] + [int(x) for x in args.steps.split(",")]
        assert steps[1] > 0

    table = runSteps(args.hal, args.outDir, args.maxBlock, args.scale,
                     steps, args.overwrite, args.maf, args.keepSequences,
                     args.trans, args.inMemory)
#    print table
    printTable(table)
if __name__ == "__main__":
    sys.exit(main())
