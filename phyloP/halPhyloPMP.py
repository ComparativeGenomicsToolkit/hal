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
import inspect
import math
import random
import string
from collections import defaultdict
from multiprocessing import Pool

from hal.stats.halStats import runParallelShellCommands
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments
from hal.stats.halStats import getHalStats
from hal.stats.halStats import getHalSequenceStats
from hal.stats.halStats import getHalGenomeLength
from hal.stats.halStats import getHalRootName


# Wrapper for hal2maf
def getHalPhyloPCmd(options):
    cmd = "halPhyloP %s %s %s %s" % (options.halFile, 
                                     options.refGenome,
                                     options.modFile,
                                     makeOutWigglePath(options))
    for opt,val in list(options.__dict__.items()):
        if (val is not None and
            (type(val) != bool or val == True) and
            (opt == 'cacheMDC' or
             opt == 'cacheRDC' or
             opt == 'cacheW0' or
             opt == 'cacheBytes' or
             opt == 'inMemory' or
             opt == 'refSequence' or
             opt == 'refTargets' or
             opt == 'start' or
             opt == 'length' or
             opt == 'targetGenomes' or
             opt == 'dupType' or
             opt == 'dupMask' or
             opt == 'step' or
             opt == 'refBed' or
             opt == 'subtree' or
             opt == 'prec')):
            if val is not True:
                cmd += ' --%s %s' % (opt, str(val))
            else:
                cmd += ' --%s' % opt
    return cmd

# Generate a WIGGLE file name from an input HAL file and some options
# the WIGGLE file will be in the specified directory
def makeOutWigglePath(options):
    wiggleFile = os.path.basename(options.wiggleFile)
    wiggleName = os.path.splitext(wiggleFile)[0]
    wiggleDir = os.path.dirname(options.wiggleFile)
    if options.sliceNumber is not None:
        wiggleName += '_%s_%d' % (options.tempID, options.sliceNumber)
    wigglePath = os.path.join(wiggleDir, wiggleName + '.wig')
    return wigglePath

# slice an interval
# return (start, length, sliceIndex)
def computeSlices(options, seqLen):
    inStart = options.start
    if inStart is None:
        inStart = 0
    inLength = options.length
    if inLength is None or inLength < 1:
        inLength = seqLen
    if inLength > 0:
        if options.sliceSize is None or options.sliceSize >= inLength:
            yield (inStart, inLength, None)
        else:
            for i in range(inLength // options.sliceSize):
                yield (inStart + i * options.sliceSize, options.sliceSize, i)
            r = inLength % options.sliceSize
            if r > 0:
                i = inLength // options.sliceSize
                yield (inStart + i * options.sliceSize, r, i)

def concatenateSlices(sliceOpts, sliceCmds):
    assert len(sliceOpts) > 0
    assert len(sliceOpts) == len(sliceCmds)
    if (sliceOpts[0].sliceSize is not None):
        for opt, cmd in zip(sliceOpts, sliceCmds):
            first = opt.sliceNumber == 0            
            sliceWigglePath = makeOutWigglePath(opt)
            assert os.path.isfile(sliceWigglePath)
            sliceNum = opt.sliceNumber
            opt.sliceNumber = None
            outWigglePath = makeOutWigglePath(opt)
            opt.sliceNumber = sliceNum
            if first:
                os.rename(sliceWigglePath, outWigglePath)
            else:
                with open(outWigglePath, "a") as tgt:
                    with open(sliceWigglePath, "r") as src:
                        for line in src:
                            tgt.write(line)
                os.remove(sliceWigglePath)

# Write out the chrom.sizes file for wigToBigWig
def writeChromSizes(options):
    if options.chromSizes is not None:
        csFile = open(options.chromSizes, "w")
        refSequenceStats = getHalSequenceStats(options.halFile, 
                                               options.refGenome)
        assert refSequenceStats is not None
        for seqStat in refSequenceStats:
            csFile.write("%s\t%s\n" % (seqStat[0], seqStat[1]))
        csFile.close()

# Decompose HAL file into slices according to the options then launch
# hal2maf in parallel processes. 
def runParallelSlices(options):
    refGenome = options.refGenome
    if refGenome is None:
        refGenome = getHalRootName(options.halFile)
    refSequenceStats = getHalSequenceStats(options.halFile, refGenome)
    options.smallFile = False
    options.firstSmallFile = True
    sliceCmds = []
    sliceOpts = []
    if options.refSequence is not None:
        refStat = [x for x in refSequenceStats if x[0] == 
                   options.refSequence]
        if len(refStat) != 1:
            raise RuntimeError("Sequence %s not found in genome %s" % (
                options.refSequence, options.refGenome))
        totalLength = int(refStat[0][1])
    else:
        totalLength = getHalGenomeLength(options.halFile, refGenome)
    
    seqOpts = copy.deepcopy(options)

    # auto compute slice size from numprocs
    if seqOpts.sliceSize == None and seqOpts.numProc > 1:
        refLen = totalLength
        if seqOpts.length is not None and seqOpts.length > 0:
            refLen = seqOpts.length
        seqOpts.sliceSize = int(math.ceil(refLen / seqOpts.numProc))
                
    for sStart, sLen, sIdx in computeSlices(seqOpts, totalLength):
        seqOpts.start = sStart
        seqOpts.length = sLen
        seqOpts.sliceNumber = sIdx
        sliceCmds.append(getHalPhyloPCmd(seqOpts))
        sliceOpts.append(copy.deepcopy(seqOpts))
            
    # run in parallel
    runParallelShellCommands(sliceCmds, options.numProc)

    # concatenate into output if desired
    concatenateSlices(sliceOpts, sliceCmds)

    writeChromSizes(options)
    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Multi-Process wrapper for halPhyloP.")

    parser.add_argument("halFile", help="Input HAL file")
    parser.add_argument("refGenome", help="Reference genome to scan")
    parser.add_argument("modFile", help="Neutral model for PhyloP. Can be "
                        "generated with halPhyloPTrain.py")
    parser.add_argument("wiggleFile", help="Output Wiggle file")

    parser.add_argument("--numProc",
                        help="Maximum number of processes to create.  If "
                        " neither --sliceSize or --splitBySequence are "
                        " specified, then the reference genome will be "
                        "sliced to require --numProc jobs",
                        type=int, default=1)
    parser.add_argument("--sliceSize",
                        help="Maximum slice of reference sequence to process "
                        "in a single process.", type=int,
                        default=None)
    parser.add_argument("--chromSizes",
                        help="Path of file to output chromosome sizes to. "
                        "Necessary for wigToBigWig",
                        default=None)

    ##################################################################
    #HDF5 OPTIONS (as copied from hal/api/hdf5_impl/hdf5CLParser.cpp)
    ##################################################################
    hdf5Grp = parser.add_argument_group('HDF5 HAL Options')
    hdf5Grp.add_argument("--cacheMDC",
                         help="number of metadata slots in hdf5 cache",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheRDC",
                         help="number of regular slots in hdf5 cache.  "
                         "should be"
                         " a prime number ~= 10 * DefaultCacheRDCBytes / "
                         "chunk",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheBytes",
                         help="maximum size in bytes of regular hdf5 cache",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheW0",
                         help="w0 parameter fro hdf5 cache", type=int,
                         default=None)
    hdf5Grp.add_argument("--inMemory",
                         help="load all data in memory (& disable hdf5 cache)",
                         action="store_true",
                         default=False)

    ##################################################################
    #HALPHYLOP OPTIONS (as copied from hal/maf/impl/hal2maf.cpp)
    ##################################################################
    hppGrp = parser.add_argument_group('halPhyloP Options')
    hppGrp.add_argument("--refSequence",
                        help="name of reference sequence within reference "
                        "genome (all sequences if empty)",
                        default=None)
    hppGrp.add_argument("--start",
                        help="coordinate within reference genome (or sequence"
                        " if specified) to start at", type=int,
                        default=None)
    hppGrp.add_argument("--length",
                        help="length of the reference genome (or sequence"
                        " if specified) to convert.  If set to 0,"
                        " the entire thing is converted", type=int,
                        default=None)
    hppGrp.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all if empty)",
                        default=None)
    hppGrp.add_argument("--dupType", 
                        help="Which duplications to mask according to dupMask "
                        "option. Choices are: "
                        "\"all\": Any duplicated region; or "
                        "\"ambiguous\": Regions within duplications where "
                        "alignments from the same species do not contain"
                        " the same base.",
                        default=None)
    hppGrp.add_argument("--dupMask",
                        help="What to do with duplicated regions. Choices are: "
                        "\"hard\": mask entire alignment column if any "
                        "duplications occur; or "
                        "\"soft\": mask species where duplications occur.",
                        default=None);
    hppGrp.add_argument("--step",
                        help="step size", type=int, default=None)
    hppGrp.add_argument("--refBed", 
                        help="Bed file with coordinates to annotate in the "
                        "reference genome to stream from standard "
                        " input.",
                        default=None)
    hppGrp.add_argument("--subtree",
                        help="Subtree root for lineage-specific acceleration/conservation",
                        default=None)
    hppGrp.add_argument("--prec",
                        help="Number of decimal places in wig output", type=int,
                        default=None)

    args = parser.parse_args()

    if not os.path.isfile(args.halFile):
        raise RuntimeError("Input hal file %s not found" % args.halFile)
    if not os.path.isfile(args.modFile):
        raise RuntimeError("Input mod file %s not found" % args.modFile)
    args.halGenomes = getHalGenomes(args.halFile)
    if not args.refGenome in args.halGenomes:
        raise RuntimeError("Reference genome %s not found." % args.refGenome)

    test = open(args.wiggleFile, "w")
    test.write("\n")
    test.close()
    os.remove(args.wiggleFile)

    if args.chromSizes is not None:
        test = open(args.chromSizes, "w")
        test.write("\n")
        test.close()
        os.remove(args.chromSizes)

    # make a little id tag for temporary slices
    S = string.ascii_uppercase + string.digits
    args.tempID = 'halPhyloPTemp' + ''.join(random.choice(S) for x in range(5))
    
    runParallelSlices(args)
    
if __name__ == "__main__":
    sys.exit(main())
