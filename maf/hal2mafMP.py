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
import inspect
import math
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
def getHal2MafCmd(options):
    cmd = "hal2maf %s %s --unique" % (options.halFile, makeOutMafPath(options))
    for opt,val in options.__dict__.items():
        if (val is not None and
            (type(val) != bool or val == True) and
            (opt == 'refGenome' or
             opt == 'refSequence' or
             opt == 'refTargets' or
             opt == 'start' or
             opt == 'length' or
             opt == 'rootGenome' or
             opt == 'targetGenomes' or
             opt == 'maxRefGap' or
             opt == 'noDupes' or
             opt == 'noAncestors' or
             opt == 'ucscNames')):
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
    halFile = os.path.basename(options.halFile)
    halName = os.path.splitext(halFile)[0]
    if options.smallFile:
        halName += '_small'
    else:
        if options.refSequence is not None:
            halName += '_%s' % options.refSequence
        if options.sliceNumber is not None:
            halName += str(options.sliceNumber)
    mafPath = os.path.join(options.mafDir, halName + '.maf')
    return mafPath

# slice an interval
# return (start, length, sliceIndex)
def computeSlices(options, seqLen):
    if seqLen > 0:
        if options.sliceSize is None or options.sliceSize >= seqLen:
            yield (0, seqLen, None)
        else:
            for i in xrange(seqLen / options.sliceSize):
                yield (i * options.sliceSize, options.sliceSize, i)
            r = seqLen % options.sliceSize
            if r > 0:
                i = seqLen / options.sliceSize
                yield (i * options.sliceSize, r, i)
    
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
    # we are going to deal with sequence coordinates
    if options.sliceBySequence is True or options.refSequence is not None:
        for sequence, seqLen, nt, nb in refSequenceStats:
            if options.refSequence is None or sequence == options.refSequence:
                seqOpts = copy.deepcopy(options)
                if seqLen < options.smallSize:
                    seqOpts.smallFile = True
                seqOpts.refGenome = refGenome
                seqOpts.refSequence = sequence
                index = 0
                for sStart, sLen, sIdx in computeSlices(seqOpts, seqLen):
                    seqOpts.start = sStart
                    seqOpts.length = sLen
                    seqOpts.sliceNumber = sIdx
                    sliceCmds.append(getHal2MafCmd(seqOpts))
                if seqOpts.smallFile is True and seqLen > 0:
                    options.firstSmallFile = False
    # we are slicing the gnome coordinates directly
    else:
        seqOpts = copy.deepcopy(options)
        genomeLen = getHalGenomeLength(seqOpts.halFile, refGenome)
        index = 0
        for sStart, sLen, sIdx in computeSlices(seqOpts, genomeLen):
            seqOpts.start = sStart
            seqOpts.length = sLen
            seqOpts.sliceNumber = sIdx
            sliceCmds.append(getHal2MafCmd(seqOpts))
            

    # run in parallel
    runParallelShellCommands(sliceCmds, options.numProc)

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Multi-Process wrapper for hal2maf.")

    parser.add_argument("halFile", help="Input HAL file")
    parser.add_argument("mafDir", help="Output directory for created MAF files")
    parser.add_argument("numProc",
                        help="Maximum number of processes to create",
                        type=int)

    parser.add_argument("--sliceSize",
                        help="Length of slice along reference genome for each "
                        "output MAF.", type=int,
                        default=None)
    parser.add_argument("--sliceBySequence",
                        help="Create an output MAF for each sequence in the "
                        "reference genome.  Note that the --sliceSize option, "
                        "if specifed, output MAF.", action="store_true",
                        default=False)
    parser.add_argument("--smallSize",
                        help="If --sliceBySequence is used, then all sequences "
                        "with length less than smallSize will be lumped into "
                        "a single output MAF called \"small\"",
                        type=int, default=0)

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
                         " a prime number ~= 10 * DefaultCacheRDCBytes / chunk",
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
    #HAL2MAF OPTIONS (as copied from hal/maf/impl/hal2maf.cpp)
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
                        help="don't write ancestral sequences",
                        action="store_true",
                        default=False)
    h2mGrp.add_argument("--ucscNames",
                        help="use UCSC convention of Genome.Seqeunce "
                        "for output names.  By default, only sequence "
                        "names are used",
                        action="store_true",
                        default=False)
    
    args = parser.parse_args()

    if not os.path.isfile(args.halFile):
        raise RuntimeError("Input hal file %s not found" % args.halFile)
    if not os.path.isdir(args.mafDir):
        os.makedirs(args.mafDir)
    if args.sliceBySequence:
        if args.start is not None:
            raise RuntimeError("--sliceBySequence option currently "
                               "incompatible with --start option")
        if args.length is not None:
            raise RuntimeError("--sliceBySequence option currently "
                               "incompatible with --length option")
    if args.sliceSize is not None and args.smallSize >= args.sliceSize:
        raise RuntimeError("--smallSize must be less than --sliceSize")

    runParallelSlices(args)
    
if __name__ == "__main__":
    sys.exit(main())
