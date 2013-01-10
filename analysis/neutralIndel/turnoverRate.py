#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Compare conserved intervals between a genome and those of its parent,
detecting patterns of conservation (overlap), gain and loss
"""
import argparse
import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np
import subprocess
import tempfile

from hal.mutations.impl.halTreeMutations import runShellCommand

def getParentGenomeName(halPath, genomeName):
    return runShellCommand("halStats %s --parent %s" % (halPath,
                                                        genomeName)).strip("\n")

# map a bed file into the parent genome's parents coordinates
def getLiftUpBedFile(halPath, genomeName, genomeBedPath, outBedPath):
    parentName = getParentGenomeName(halPath, genomeName)
    runShellCommand("halLiftover %s %s %s %s %s" % (halPath,
                                                    genomeName,
                                                    genomeBedPath,
                                                    parentName,
                                                    outBedPath))
                    
# compute aligned regions of a genome as a bed
def getAlignedBed(halPath, genomeName, outBedPath):
    runShellCommand("halAlignedExtract %s %s --alignedFile %s" % (halPath,
                                                                  genomeName,
                                                                  outBedPath))
    
# compute inbed1 AND inbed1                    
def getIntersectBed(inBed1, inBed2, outBed):
    runShellCommand("intersectBed -a %s -b %s > %s" % (inBed1, inBed2, outBed))

# compute inbed1 - inbed2
def getSubtractBed(inBed1, inBed2, outBed):
    runShellCommand("subtractBed -a %s -b %s > %s" % (inBed1, inBed2, outBed))

# compute inbed1 OR inbed2
def getUnionBed(inBed1, inBed2, outBed):
    runShellCommand("cat %s %s | sortBed | mergeBed > %s" % (inBed1, inBed2,
                                                             outBed))

# compute total length of all intervals in a bed file.  NOTE THAT
# OVERLAPPING OR DUPLICATE INTERVALS WILL BE COUNTED AS WELL, SO MAKE SURE
# INPUT DOES NOT CONTAIN OVERLAPS!
def getBedLength(bedPath):
    length = 0
    bedFile = open(bedPath)
    for line in bedFile:
        clnLine = line.strip()
        if len(clnLine) > 0 and clnLine[0] != "#":
            toks = clnLine.split()
            if len(toks) > 2:
                lineLength = int(toks[2]) - int(toks[1])
                length += lineLength
    bedFile.close()
    return length
        
# do analysis along a branch.  Given conservaiton beds for the genome and
# its parent, do a basic set breakdown to find the intersection and differences
# which correspond respectively to conservation, and gain/loss.
def compareConservationOverBranch(halPath, genomeName, genomeBed, parentBed,
                                  outMappedGenomeBed, outConservationBed,
                                  outAlignedBed,
                                  outGainBed, outLossBed):
    #1) map genome to parent's coordinates
    getLiftUpBedFile(halPath, genomeName, genomeBed, outMappedGenomeBed)
    
    #2) compute intersection
    getIntersectBed(parentBed, outMappedGenomeBed, outConservationBed)

    #3) compute gain
    getSubtractBed(outMappedGenomeBed, parentBed, outGainBed)
    
    #4) compute loss
    getSubtractBed(parentBed, outMappedGenomeBed, outLossBed)

    #5) compute aligned regions                
    getAlignedBed(halPath, genomeName, outAlignedBed)
    
    conLen = getBedLength(outConservationBed)
    gainLen = getBedLength(outGainBed)
    lossLen = getBedLength(outLossBed)
    unconLen = getBedLength(outAlignedBed) - conLen - gainLen - lossLen

    return (conLen, gainLen, lossLen, unconLen)
                                                    
                                                    
            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("halPath", type=str,
                        help="path to hal file")
    parser.add_argument("genome", type=str,
                        help="name of genome")
    parser.add_argument("--bed", type=str,
                        help="bed file for conserved regions "
                        "[default=genome_cons.bed]",
                        default=None)
    parser.add_argument("--parBed", type=str,
                        help="bed file conserved regions in parent "
                        "[default=parentGenome_cons.bed]",
                        default=None)
    parser.add_argument("--workDir", type=str,
                        help="default dir for bed files", default=None)
    
    args = parser.parse_args()
    workDir = args.workDir
    if workDir is None:
        if args.parBed is not None:
            workDir = os.path.dirname(args.parBed)
        elif args.bed is not None:
            workDir = os.path.dirname(args.bed)
        else:
            raise RuntimeError("need to specify at least one of "
                               "--bed, --parBed, or --workDir")
    
    genomeBed = args.bed
    if genomeBed is None:
        genomeBed = os.path.join(workDir, args.genome + "_cons.bed")
    parentBed = args.parBed
    if parentBed is None:
        parentBed = os.path.join(workDir, getParentGenomeName(args.halPath,
                                                              args.genome) +
                                 "_cons.bed")
    outMappedGenomeBed = os.path.join(workDir, args.genome + "_pm.bed")
    outConservationBed = os.path.join(workDir, args.genome + "_int.bed")
    outAlignedBed = os.path.join(workDir, args.genome + "_al.bed")
    outGainBed = os.path.join(workDir, args.genome + "_gain.bed")
    outLossBed = os.path.join(workDir, args.genome + "_loss.bed")

    (conLen, gainLen, lossLen, unconLen) = compareConservationOverBranch(
        args.halPath, args.genome, genomeBed, parentBed,
        outMappedGenomeBed, outConservationBed, outAlignedBed,
        outGainBed, outLossBed)

    print (conLen, gainLen, lossLen, unconLen)
    
if __name__ == "__main__":
    sys.exit(main())

