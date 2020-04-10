#!/usr/bin/env python3

# Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
# Released under the MIT license, see LICENSE.txt

"""Run hal2maf in parallel by slicing along reference genome.
"""
import argparse
import os
import sys
import copy
import random
import math
import statistics
import string
from glob import glob
from sonLib.bioio import system

from hal.stats.halStats import runParallelShellCommands
from hal.stats.halStats import getHalSequenceStats
from hal.stats.halStats import getHalGenomeLength
from hal.stats.halStats import getHalRootName


# Wrapper for hal2maf
def getHal2MafCmd(options):
    cmd = "hal2maf %s %s --unique" % (options.halFile, makeOutMafPath(options))
    for opt, val in list(options.__dict__.items()):
        if (val is not None and
            ((not isinstance(val, bool) or (val == True)) and
             (opt in ('cacheMDC', 'cacheRDC', 'cacheW0', 'cacheBytes', 'hdf5InMemory', 'refGenome',
                      'refSequence', 'refTargets', 'start', 'length', 'rootGenome',
                      'targetGenomes', 'maxRefGap', 'noDupes', 'noAncestors', 'onlySequenceNames')))):
            if val is not True:
                cmd += ' --%s %s' % (opt, str(val))
            else:
                cmd += ' --%s' % opt
    if options.smallFile and not options.firstSmallFile:
        cmd += ' --append'
    return cmd

# Generate a MAF file name from an input HAL file and some options
# the MAF file will be in the specified directory
def makeOutMafPath(options):
    mafFile = os.path.basename(options.mafFile)
    mafName = os.path.splitext(mafFile)[0]

    mafDir = os.path.dirname(options.mafFile)
    if options.smallFile:
        mafName += '_small'
    else:
        if options.refSequence is not None and options.splitBySequence is True:
            mafName += '_%s' % options.refSequence
        if options.sliceNumber is not None:
            mafName += '_%s_%d' % (options.tempID, options.sliceNumber)
    mafPath = os.path.join(mafDir, mafName + '.maf')
    return mafPath

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
    for opt, cmd in zip(sliceOpts, sliceCmds):
        first = opt.sliceNumber == 0
        sliceMafPath = makeOutMafPath(opt)
        if os.path.isfile(sliceMafPath) and opt.sliceNumber is not None:
            sliceNum = opt.sliceNumber
            opt.sliceNumber = None
            outMafPath = makeOutMafPath(opt)
            opt.sliceNumber = sliceNum
            if first:
                os.rename(sliceMafPath, outMafPath)
            else:
                with open(outMafPath, "a") as tgt:
                    with open(sliceMafPath, "r") as src:
                        for line in src:
                            if not line[0] == '#':
                                tgt.write(line)
                os.remove(sliceMafPath)

def splitBed(bed, numParts):
    """Split up a bed file by lines into N parts, return the paths of the
    split files"""
    with open(bed) as f:
        numLines = sum(1 for line in f)
    # Random suffix so two runs on the same file don't collide
    suffix = "".join([random.choice(string.ascii_uppercase) for _ in range(7)])
    system("split -l %d %s %s.temp.%s" % (math.ceil(float(numLines)/numParts),
                                          bed, bed, suffix))
    return glob('%s.temp.%s*' % (bed, suffix))

def partitionRefTargets(options):
    "partition by refTargets option"
    splitBedFiles = splitBed(options.refTargets, options.numProc)
    sliceCmds = []
    sliceOpts = []
    for i, splitBedFile in enumerate(splitBedFiles):
        seqOpts = copy.deepcopy(options)
        seqOpts.refTargets = splitBedFile
        seqOpts.sliceNumber = i
        sliceCmds.append(getHal2MafCmd(seqOpts))
        sliceOpts.append(seqOpts)
    return sliceCmds, sliceOpts

def partitionBySeqCoords(options, refGenome):
    "we are going to deal with sequence coordinates"
    refSequenceStats = getHalSequenceStats(options.halFile, refGenome)
    sliceCmds = []
    sliceOpts = []
    for sequence, seqLen, nt, nb in refSequenceStats:
        if options.refSequence is None or sequence == options.refSequence:
            seqOpts = copy.deepcopy(options)
            if seqLen < options.smallSize:
                seqOpts.smallFile = True
            seqOpts.refGenome = refGenome
            seqOpts.refSequence = sequence
            for sStart, sLen, sIdx in computeSlices(seqOpts, seqLen):
                seqOpts.start = sStart
                seqOpts.length = sLen
                seqOpts.sliceNumber = sIdx
                sliceCmds.append(getHal2MafCmd(seqOpts))
                sliceOpts.append(copy.deepcopy(seqOpts))
            if seqOpts.smallFile is True and seqLen > 0:
                options.firstSmallFile = False
    return sliceCmds, sliceOpts

def partitionByGenomeCoords(options, refGenome):
    "we are slicing the gnome coordinates directly"
    sliceCmds = []
    sliceOpts = []
    seqOpts = copy.deepcopy(options)
    assert seqOpts.splitBySequence is False
    refSequenceStats = getHalSequenceStats(options.halFile, refGenome)
    # auto compute slice size from numprocs
    if seqOpts.sliceSize is None and seqOpts.numProc > 1:
        if seqOpts.length is not None and seqOpts.length > 0:
            refLen = seqOpts.length
        else:
            # use median of sequence lengths
            refLen = int(statistics.median([r[1] for r in refSequenceStats]))
        seqOpts.sliceSize = math.ceil(math.ceil(refLen / seqOpts.numProc))

    for refSeqStat in refSequenceStats:
        seqOpts.refSequence = refSeqStat[0]
        for sStart, sLen, sIdx in computeSlices(seqOpts, refSeqStat[1]):
            seqOpts.start = sStart
            seqOpts.length = sLen
            seqOpts.sliceNumber = sIdx
            sliceCmds.append(getHal2MafCmd(seqOpts))
            sliceOpts.append(copy.deepcopy(seqOpts))
    return sliceCmds, sliceOpts

# Decompose HAL file into slices according to the options then launch
# hal2maf in parallel processes.
def runParallelSlices(options):
    refGenome = options.refGenome
    if refGenome is None:
        refGenome = getHalRootName(options.halFile)
    options.smallFile = False
    options.firstSmallFile = True
    if options.refTargets:
        sliceCmds, sliceOpts = partitionRefTargets(options)
    elif options.splitBySequence is True or options.refSequence is not None:
        sliceCmds, sliceOpts = partitionBySeqCoords(options, refGenome)
    else:
        sliceCmds, sliceOpts = partitionByGenomeCoords(options, refGenome)

    # run in parallel
    runParallelShellCommands(sliceCmds, options.numProc)

    # clean up temporary bed files (if present)
    for opts in sliceOpts:
        if opts.refTargets and os.path.isfile(opts.refTargets):
            os.remove(opts.refTargets)

    # concatenate into output if desired
    concatenateSlices(sliceOpts, sliceCmds)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Multi-Process wrapper for hal2maf.")

    parser.add_argument("halFile", help="Input HAL file")
    parser.add_argument("mafFile", help="Output MAF file")

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
    parser.add_argument("--splitBySequence",
                        help="Create an output MAF for each sequence in the "
                        "reference genome. Output files will have the format "
                        "mafFile_sequenceName.maf", action="store_true",
                        default=False)
    parser.add_argument("--smallSize",
                        help="If --splitBySequence is used, then all sequences "
                        "with length less than smallSize will be lumped into "
                        "a single output MAF called \"mafFile_small.maf\"",
                        type=int, default=0)

    ##################################################################
    # HDF5 OPTIONS (as copied from hal/api/hdf5_impl/hdf5CLParser.cpp)
    ##################################################################
    hdf5Grp = parser.add_argument_group('HDF5 HAL Options')
    hdf5Grp.add_argument("--cacheMDC",
                         help="number of metadata slots in hdf5 cache",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheRDC",
                         help="number of regular slots in hdf5 cache.  "
                         "should be"
                         " a prime number ~= 10 * DefaultCacheRDCBytes / chunk",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheBytes",
                         help="maximum size in bytes of regular hdf5 cache",
                         type=int,
                         default=None)
    hdf5Grp.add_argument("--cacheW0",
                         help="w0 parameter fro hdf5 cache", type=float,
                         default=None)
    hdf5Grp.add_argument("--hdf5InMemory",
                         help="load all data in memory (& disable hdf5 cache)",
                         action="store_true",
                         default=False)

    ##################################################################
    # HAL2MAF OPTIONS (as copied from hal/maf/impl/hal2maf.cpp)
    ##################################################################
    h2mGrp = parser.add_argument_group('hal2maf Options')
    h2mGrp.add_argument("--refGenome",
                        help="name of reference genome (root if empty)",
                        default=None)
    h2mGrp.add_argument("--refSequence",
                        help="name of reference sequence within reference "
                        "genome (all sequences if empty)",
                        default=None)
    h2mGrp.add_argument("--refTargets",
                        help="bed file coordinates of intervals in the"
                        " reference genome to export",
                        default=None)
    h2mGrp.add_argument("--start",
                        help="coordinate within reference genome (or sequence"
                        " if specified) to start at", type=int,
                        default=None)
    h2mGrp.add_argument("--length",
                        help="length of the reference genome (or sequence"
                        " if specified) to convert.  If set to 0,"
                        " the entire thing is converted", type=int,
                        default=None)
    h2mGrp.add_argument("--rootGenome",
                        help="name of root genome (none if empty)",
                        default=None)
    h2mGrp.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all if empty)",
                        default=None)
    h2mGrp.add_argument("--maxRefGap",
                        help="maximum gap length in reference", type=int,
                        default=None)
    h2mGrp.add_argument("--noDupes",
                        help="ignore paralogy edges",
                        action="store_true",
                        default=False)
    h2mGrp.add_argument("--noAncestors",
                        help="don't write ancestral sequences. IMPORTANT: "
                        "Must be used in conjunction with --refGenome"
                        " to set a non-ancestral genome as the reference"
                        " because the default reference is the root.",
                        action="store_true",
                        default=False)
    h2mGrp.add_argument("--onlySequenceNames",
                        help="use sequence names "
                        "for output names.  By default, the UCSC "
                        "convention of Genome.Sequence is used",
                        action="store_true",
                        default=False)

    args = parser.parse_args()
    if args.refTargets and any([args.splitBySequence, args.sliceSize, args.start, args.length, args.refSequence]):
        parser.error("--refTargets not compatible with --splitBySequence, --sliceSize, --start, --length, or --refSequence")
    if args.splitBySequence:
        if args.start is not None:
            parser.error("--splitBySequence option currently incompatible with --start option")
        if args.length is not None:
            parser.error("--splitBySequence option currently incompatible with --length option")
    if args.sliceSize is not None and args.smallSize >= args.sliceSize:
        parser.error("--smallSize must be less than --sliceSize")

    if not os.path.isfile(args.halFile):
        raise RuntimeError("Input hal file %s not found" % args.halFile)
    assert args.mafFile != args.halFile
    test = open(args.mafFile, "w")
    test.write("\n")
    test.close()
    os.remove(args.mafFile)
    # make a little id tag for temporary maf slices
    S = string.ascii_uppercase + string.digits
    args.tempID = 'hal2mafTemp' + ''.join(random.choice(S) for x in range(5))
    runParallelSlices(args)

if __name__ == "__main__":
    sys.exit(main())
