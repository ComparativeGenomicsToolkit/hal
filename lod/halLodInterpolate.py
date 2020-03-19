#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Generate a series of HAL files at progressively coarse levels of detail
from an input file by calling halLodExtract
"""
import argparse
import os
import sys
import copy
import subprocess
import time
import math
from collections import defaultdict
from multiprocessing import Pool

from hal.stats.halStats import runParallelShellCommands
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments
from hal.stats.halStats import getHalStats
from hal.stats.halStats import getHalSequenceStats

# specify upper limit of lods.
# (MUST MANUALLY KEEP CONSISTENT WITH global LodManager::MaxLodToken
# variable in hal/lod/impl/halLodManager.cpp)
MaxLodToken = "max"

# Wrapper for halLodExtract
def getHalLodExtractCmd(inHalPath, outHalPath, scale, keepSeq, inMemory,
                     probeFrac, minSeqFrac, chunk, minCovFrac):
    cmd = "halLodExtract %s %s %s" % (inHalPath, outHalPath, scale)
    if keepSeq is True:
        cmd += " --keepSequences"
    if inMemory is True:
        cmd += " --inMemory"
    if probeFrac is not None:
        cmd += " --probeFrac %f" % probeFrac
    if minSeqFrac is not None:
        cmd += " --minSeqFrac %f" % minSeqFrac
    if chunk is not None and chunk > 0:
        cmd += " --chunk %d" % chunk

    return cmd

# All created paths get put in the same place using the same logic
def makePath(inHalPath, outDir, step, name, ext):
    inFileName = os.path.splitext(os.path.basename(inHalPath))[0]
    outPath = os.path.join(outDir, "%s_%s_%i.%s" % (inFileName,
                                                    name, step, ext))
    return outPath

# Return the length of the longest genome in the HAL file
def getMaxGenomeLength(statsTable):
    maxLen = 0
    for row in statsTable:
        maxLen = max(maxLen, int(row[2]))
    return maxLen

# Return the smallest averege block length of any genome in the HAL file
def getMinAvgBlockSize(statsTable):
    minAvgBlockSize = sys.maxsize
    for row in statsTable:
        if float(row[3]) > 0:
            avgTop = float(row[2]) / float(row[3])
            minAvgBlockSize = min(minAvgBlockSize, avgTop)
        if float(row[4]) > 0:
            avgBottom = float(row[2]) / float(row[4])
            minAvgBlockSize = min(minAvgBlockSize, avgBottom)
    assert minAvgBlockSize > 0 and minAvgBlockSize != sys.maxsize
    return minAvgBlockSize

# Return the smallest fraction of any genome that would be left if we
# cut out all sequences of length less than minLen
def getMinCoverageFrac(sequenceStatsTable, minLen):
    minCoverage = 1.0
    for genome, sequenceStats in list(sequenceStatsTable.items()):
        totalLength = 0.0
        uncutLength = 0.0
        for sequence, seqLen, numTop, numBot in sequenceStats:
            totalLength += seqLen
            if seqLen >= minLen:
                uncutLength += seqLen
        coverage = uncutLength / totalLength
        minCoverage = min(coverage, minCoverage)

    return minCoverage        
        
# Get a lest of step-sizes required to interpolate the hal such that
# the maximum level of detail has at most maxBlock (expected) blocks
# per genome, and each hal file has multFac (approx maxBlock) bigger
# blocks than the previous
def getSteps(halPath, maxBlock, scaleFactor, minLod0, cutOffFrac, minSeqFrac,
            minCovFrac):
    statsTable = getHalStats(halPath)
    sequenceStatsTable = dict()
    for row in statsTable:
        sequenceStatsTable[row[0]] = getHalSequenceStats(halPath, row[0])
    maxLen = getMaxGenomeLength(statsTable)
    assert maxLen > 0
    maxStep = math.ceil(float(maxLen) / float(maxBlock))
    lodBaseStep = math.ceil(float(minLod0) / float(maxBlock))
    baseStep = max(lodBaseStep, getMinAvgBlockSize(statsTable))
    outList = []
    step = baseStep
    # last LOD is just "max" token which tells browser it and anything
    # beyond is disabled.
    lastIsMax = False
    while True:
        outList.append(step)
        if step > maxStep * cutOffFrac:
            break
        minCoverage = 1.0
        if minSeqFrac > 0. and minCovFrac > 0.:
            minCoverageFrac = getMinCoverageFrac(sequenceStatsTable,
                                                 math.floor(step * minSeqFrac))
            if minCoverageFrac < minCovFrac:
                lastIsMax = True
                break
        step *= scaleFactor
    return [int(x) for x in outList], lastIsMax

def formatOutHalPath(outLodPath, outHalPath, absPath):
    if absPath:
        return os.path.abspath(outHalPath)
    else:
        return os.path.relpath(outHalPath, os.path.dirname(outLodPath))
                         
# Run halLodExtract for each level of detail.
def createLods(halPath, outLodPath, outDir, maxBlock, scale, overwrite,
               maxDNA, absPath, trans, inMemory, probeFrac, minSeqFrac,
               scaleCorFac, numProc, chunk, minLod0, cutOff, minCovFrac):
    lodFile = open(outLodPath, "w")
    lodFile.write("0 %s\n" % formatOutHalPath(outLodPath, halPath, absPath))
    steps, lastIsMax = getSteps(halPath, maxBlock, scale, minLod0, cutOff,
                                minSeqFrac, minCovFrac)
    curStepFactor = scaleCorFac
    lodExtractCmds = []
    prevStep = None
    for stepIdx in range(1,len(steps)):
        step = int(max(1, steps[stepIdx] * curStepFactor))
        maxQueryLength = maxBlock * steps[stepIdx - 1]
        keepSequences = maxQueryLength <= maxDNA
        #we no longer pass the step to the halLodExtract executable,
        #rather we give the corresponding the scale factor and let
        #the step get computed for each internal node (instead of using the step
        #here which is a global minimum
        stepScale = (scale ** stepIdx) * curStepFactor
        outHalPath = makePath(halPath, outDir, step, "lod", "hal")
        srcPath = halPath
        if trans is True and stepIdx > 1:
            srcPath = makePath(halPath, outDir, prevStep, "lod", "hal")
        isMaxLod = stepIdx == len(steps) - 1 and lastIsMax is True
        if not isMaxLod and (overwrite is True or
                             not os.path.isfile(outHalPath)):
            lodExtractCmds.append(
                getHalLodExtractCmd(srcPath, outHalPath, stepScale,
                                    keepSequences, inMemory, probeFrac,
                                    minSeqFrac, chunk, minCovFrac))
        lodPath =  formatOutHalPath(outLodPath, outHalPath, absPath)
        if isMaxLod:
            lodPath = MaxLodToken
        
        lodFile.write("%d %s\n" % (maxQueryLength, lodPath))

        if prevStep is not None and prevStep > steps[-1]:
            break
        prevStep = step
        curStepFactor *= scaleCorFac
    lodFile.close()
    runParallelShellCommands(lodExtractCmds, numProc)
    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate a series of HAL files at progressively coarse "
        "levels of detail from an input file by calling halLodExtract")

    parser.add_argument("hal", help="input hal")
    parser.add_argument("outLodFile", help="output text file with links to" 
                        " interpolated hal files.  with each file is"
                        " associated a value stating its minimum "
                        "suggested query range (in bases)")
    parser.add_argument("--maxBlock",
                        help="maximum desired number of blocks to ever " 
                        "display at once.", type=int,
                        default=223)
    parser.add_argument("--scale",
                        help="scaling factor between two successive levels"
                        " of detail", type=float,
                        default=2.5)
    parser.add_argument("--outHalDir", help="path of directory where "
                        "interpolated hal files are stored.  By default "
                        "they will be stored in the same directory as the "
                        "input file",
                        default=None)
    parser.add_argument("--resume",
                        help="do not overwrite existing hal lod files if they "
                        "exist.",
                        action="store_true", default=False)
    parser.add_argument("--maxDNA",
                        help="maximum DNA sequence query.  Generated levels of"
                        " detail with associated minimum query ranges > maxDNA"
                        " will not contain sequence information.  -1 can be "
                        "used to specify that all levels will get sequence",
                        type=int,
                        default=0)
    parser.add_argument("--absPath",
                        help="write absolute path of created HAL files in the"
                        " outLodFile.  By default, the paths are relative to "
                        "the path of the outLodFile.",
                        action="store_true", default=False)
    parser.add_argument("--trans", help="Generate level of detail X from "
                        "X-1.  By default, all levels of detail are generated "
                        "from the original HAL (X=0)",
                        action="store_true", default=False)
    parser.add_argument("--inMemory", help="Load entire hdf5 arrays into "
                        "memory, overriding cache.",
                        action="store_true", default=False)
    parser.add_argument("--probeFrac", help="Fraction of bases in step-interval"
                        " to sample while looking for most aligned column. "
                        "Use default from halLodExtract if not set.  To see"
                        " default value, use halLodExtract --help",                
                        type=float, default=None)
    parser.add_argument("--minSeqFrac", help="Minumum sequence length to sample "
                        "as fraction of step size: ie sequences with "
                        "length <= floor(minSeqFrac * step) are ignored."
                        "Use default from halLodExtract if not set. To see"
                        " default value, use halLodExtract --help",                  
                        # Note: needs to be manually synched with 
                        # value in halLodInterpolate.py
                        type=float, default=0.5)
    parser.add_argument("--minCovFrac", help="Minimum fraction of a genome"
                        " that must be covered by sequences that exceed"
                        " --minSeqFrac * step.  LODs that would violate this"
                        " threshold will not be generated (or displayed in"
                        " the browser).  This is seen a better than the "
                        "alternative, which is to produce unreasonably sparse"
                        " LODs because half the sequences were not sampled",
                        type=float, default=0.9)
    parser.add_argument("--scaleCorFac", help="Correction factor for scaling. "
                        " Assume that scaling by (X * scaleCorFactor) is "
                        " required to reduce the number of blocks by X.",
                        type=float, default=1.0)
    parser.add_argument("--numProc", help="Number of concurrent processes",
                        type=int, default=1)
    parser.add_argument("--chunk", help="Chunk size of output hal files.  ",
                        type=int, default=None)
    parser.add_argument("--minLod0", help="Override other parameters to "
                        "ensure that Lod0 (original hal) gets range from 0 "
                        "to at least the specified value", type=int,
                        default=0)
    parser.add_argument("--cutOff", help="Cut-off fraction for determining "
                        "highest level of detail.  Normally, --maxBlocks is"
                        " used to determine the step-size for each LOD, but "
                        "the exponential scaling can cause this to create"
                        " final (highest) LODs that are too sparse.  The lower"
                        "this parameter is (must stay > 0), the sooner "
                        "the LOD generation process will be cut off, and the"
                        " more fine-grained the highest LOD will be",
                        default=0.75, type=float)

    args = parser.parse_args()

    if args.outHalDir is not None and not os.path.exists(args.outHalDir):
        os.makedirs(args.outHalDir)

    if not os.path.isfile(args.hal):
        raise RuntimeError("Input hal file %s not found" % args.hal)
    if args.outHalDir is None:
        args.outHalDir = os.path.dirname(os.path.abspath(args.hal))
    if not os.path.isdir(args.outHalDir):
        raise RuntimeError("Invalid output directory %s" % args.outHalDir)
    assert args.scaleCorFac > 0
    if args.trans is True and args.numProc > 1:
        raise RuntimeError("--numProc > 1 not supported when --trans option is "
                           "set")

    if args.maxDNA < 0:
        args.maxDNA = sys.maxsize

    createLods(args.hal, args.outLodFile, args.outHalDir,
               args.maxBlock, args.scale, not args.resume, args.maxDNA,
               args.absPath, args.trans, args.inMemory, args.probeFrac,
               args.minSeqFrac, args.scaleCorFac, args.numProc, args.chunk,
               args.minLod0, args.cutOff, args.minCovFrac)
    
if __name__ == "__main__":
    sys.exit(main())
