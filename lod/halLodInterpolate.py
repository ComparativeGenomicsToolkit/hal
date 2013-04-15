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

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments
from hal.stats.halStats import getHalStats

from hal.lod.halLodBenchmark import runHalLodExtract
from hal.lod.halLodBenchmark import makePath

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

# Try to find a factor smaller than inFactor such that the range is evenly
# divided into exponents of the factor
def getFactor(first, last, inFactor):
    assert first >= 0 and last > first
    assert inFactor > 0
    distance = (last - first) + 1
    tgtRes = math.ceil(math.log(distance, inFactor))
    bestFactor = None
    bestDelta = sys.maxint
    for testFactor in xrange(inFactor, 2, -1):
        res = math.log(distance, testFactor)
        if res <= tgtRes and tgtRes - res < bestDelta:
            bestFactor = testFactor
            bestDelta = tgtRes -res
    assert bestFactor > 1
    return bestFactor
        
# Get a lest of step-sizes required to interpolate the hal such that
# the maximum level of detail has at most maxBlock (expected) blocks
# per genome, and each hal file has multFac (approx maxBlock) bigger
# blocks than the previous
def getSteps(halPath, maxBlock):
    maxLen = getMaxGenomeLength(halPath)
    assert maxLen > 0
    maxStep = int(math.ceil(float(maxLen) / float(maxBlock)))
    baseStep = int(getMinAvgBlockSize(halPath))
    multFac = getFactor(baseStep, maxStep, maxBlock)
    outList = []
    step = baseStep * multFac
    while step * multFac < maxLen:
        outList.append((step, step*multFac))
        step *= multFac
    return outList

# Run halLodExtract for each level of detail.
def createLods(halPath, outLodPath, outDir, maxBlock, overwrite,
               keepSequences):
    lodFile = open(outLodPath, "w")
    lodFile.write("0 %s\n" % halPath)
    for iteration in getSteps(halPath, maxBlock):
        step = iteration[0]
        blen = iteration[1]
        outHalPath = makePath(halPath, outDir, step, "lod", "hal")
        if overwrite is True or not os.path.isfile(outHalPath):
            runHalLodExtract(halPath, outHalPath, step, keepSequences)
        lodFile.write("%d %s\n" % (blen, outHalPath))
    lodFile.close()
    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Generate a series of HAL files at progressively coarse "
        "levels of detail from an input file by calling halLodExtract")

    parser.add_argument("hal", help="input hal")
    parser.add_argument("outLodFile", help="output text file with links to" 
                        " interpolated hal files.")
    parser.add_argument("--maxBlock",
                        help="maximum desired number of blocks to ever " 
                        "display at once.", type=int,
                        default=1000)
    parser.add_argument("--outHalDir", help="path of directory where "
                        "interpolated hal files are stored.  By default "
                        "they will be stored in the same directory as the "
                        "input file",
                        default=None)
    parser.add_argument("--overwrite",
                        help="overwrite existing hal files if they exist.",
                        action="store_true", default=False)
    parser.add_argument("--keepSequences",
                        help="copy DNA sequences into interpolated HAL files",
                        action="store_true", default=False)
        
    args = parser.parse_args()

    if args.outHalDir is not None and not os.path.exists(args.outHalDir):
        os.makedirs(args.outHalDir)

    if not os.path.isfile(args.hal):
        raise RuntimeError("Input hal file %s not found" % args.hal)
    if args.outHalDir is None:
        args.outHalDir = os.path.dirname(args.hal)
    if not os.path.isdir(args.outHalDir):
        raise RuntimeError("Invalid output directory %s" % args.outHalDir)

    createLods(args.hal, args.outLodFile, args.outHalDir,
               args.maxBlock, args.overwrite, args.keepSequences)
    
if __name__ == "__main__":
    sys.exit(main())
