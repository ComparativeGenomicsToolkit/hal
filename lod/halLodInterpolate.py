#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

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

# Wrapper for halLodExtract
def getHalLodExtractCmd(inHalPath, outHalPath, step, keepSeq, inMemory,
                     probeFrac, minSeqFrac, chunk):
    cmd = "halLodExtract %s %s %s" % (inHalPath, outHalPath, step)
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
def getMaxGenomeLength(halPath):
    statsTable = getHalStats(halPath)
    maxLen = 0
    for row in statsTable:
        maxLen = max(maxLen, int(row[2]))
    return maxLen

# Return the smallest averege block length of any genome in the HAL file
def getMinAvgBlockSize(halPath):
    statsTable = getHalStats(halPath)
    minAvgBlockSize = sys.maxint
    for row in statsTable:
        if float(row[3]) > 0:
            avgTop = float(row[2]) / float(row[3])
            minAvgBlockSize = min(minAvgBlockSize, avgTop)
        if float(row[4]) > 0:
            avgBottom = float(row[2]) / float(row[4])
            minAvgBlockSize = min(minAvgBlockSize, avgBottom)
    assert minAvgBlockSize > 0 and minAvgBlockSize != sys.maxint
    return minAvgBlockSize
        
# Get a lest of step-sizes required to interpolate the hal such that
# the maximum level of detail has at most maxBlock (expected) blocks
# per genome, and each hal file has multFac (approx maxBlock) bigger
# blocks than the previous
def getSteps(halPath, maxBlock, scaleFactor):
    maxLen = getMaxGenomeLength(halPath)
    assert maxLen > 0
    maxStep = math.ceil(float(maxLen) / float(maxBlock))
    baseStep = getMinAvgBlockSize(halPath)
    outList = []
    step = baseStep
    while True:
        outList.append(step)
        if step > maxStep:
            break
        step *= scaleFactor
    return [int(x) for x in outList]

def formatOutHalPath(outLodPath, outHalPath, absPath):
    if absPath:
        return os.path.abspath(outHalPath)
    else:
        return os.path.relpath(outHalPath, os.path.dirname(outLodPath))
                         
# Run halLodExtract for each level of detail.
def createLods(halPath, outLodPath, outDir, maxBlock, scale, overwrite,
               maxDNA, absPath, trans, inMemory, probeFrac, minSeqFrac,
               scaleCorFac, numProc, chunk):
    lodFile = open(outLodPath, "w")
    lodFile.write("0 %s\n" % formatOutHalPath(outLodPath, halPath, absPath))
    steps = getSteps(halPath, maxBlock, scale)
    curStepFactor = scaleCorFac
    lodExtractCmds = []
    prevStep = None
    for stepIdx in xrange(1,len(steps)):
        step = int(max(1, steps[stepIdx] * curStepFactor))
        maxQueryLength = maxBlock * steps[stepIdx - 1]
        keepSequences = maxQueryLength <= maxDNA
        outHalPath = makePath(halPath, outDir, step, "lod", "hal")
        srcPath = halPath
        if trans is True and stepIdx > 1:
            srcPath = makePath(halPath, outDir, prevStep, "lod", "hal")  
        if overwrite is True or not os.path.isfile(outHalPath):
            lodExtractCmds.append(
                getHalLodExtractCmd(srcPath, outHalPath, step, keepSequences,
                                    inMemory, probeFrac, minSeqFrac, chunk))
  
        lodFile.write("%d %s\n" % (maxQueryLength,
                                   formatOutHalPath(outLodPath, outHalPath,
                                                    absPath)))
        if prevStep > steps[-1]:
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
                        default=123)
    parser.add_argument("--scale",
                        help="scaling factor between two successive levels"
                        " of detail", type=float,
                        default=3.0)
    parser.add_argument("--outHalDir", help="path of directory where "
                        "interpolated hal files are stored.  By default "
                        "they will be stored in the same directory as the "
                        "input file",
                        default=None)
    parser.add_argument("--overwrite",
                        help="overwrite existing hal files if they exist.",
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
                        "Use default from halLodExtract if not set.",                
                        type=float, default=None)
    parser.add_argument("--minSeqFrac", help="Minumum sequence length to sample "
                        "as fraction of step size: ie sequences with "
                        "length <= floor(minSeqFrac * step) are ignored."
                        "Use default from halLodExtract if not set.",                
                        type=float, default=None)
    parser.add_argument("--scaleCorFac", help="Correction factor for scaling. "
                        " Assume that scaling by (X * scaleCorFactor) is "
                        " required to reduce the number of blocks by X.",
                        type=float, default=1.3)
    parser.add_argument("--numProc", help="Number of concurrent processes",
                        type=int, default=1)
    parser.add_argument("--chunk", help="Chunk size of output hal files.  "
                        "Note we keep the default lower than regular HAL "
                        "(1000) to try to favor network access",
                        type=int, default=200)

    args = parser.parse_args()

    if args.outHalDir is not None and not os.path.exists(args.outHalDir):
        os.makedirs(args.outHalDir)

    if not os.path.isfile(args.hal):
        raise RuntimeError("Input hal file %s not found" % args.hal)
    if args.outHalDir is None:
        args.outHalDir = os.path.dirname(args.hal)
    if not os.path.isdir(args.outHalDir):
        raise RuntimeError("Invalid output directory %s" % args.outHalDir)
    assert args.scaleCorFac > 0
    if args.trans is True and args.numProc > 1:
        raise RuntimeError("--numProc > 1 not supported when --trans option is "
                           "set")

    if args.maxDNA < 0:
        args.maxDNA = sys.maxint

    createLods(args.hal, args.outLodFile, args.outHalDir,
               args.maxBlock, args.scale, args.overwrite, args.maxDNA,
               args.absPath, args.trans, args.inMemory, args.probeFrac,
               args.minSeqFrac, args.scaleCorFac, args.numProc, args.chunk)
    
if __name__ == "__main__":
    sys.exit(main())
