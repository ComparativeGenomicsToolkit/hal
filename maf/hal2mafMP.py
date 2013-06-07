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

# Wrapper for hal2maf
def getHal2MafCmd(options):
    cmd = "hal2maf %s %s" % (options.inHalPath, options.outHalPath)
    for opt,val in options:
        if (val is not None and
            opt != 'halFile' and
            opt != 'mafDir' and
            opt != 'numProc' and
            opt != 'sliceSize' and
            opt != 'sliceBySequence'):
            cmd += ' %s %s' % (opt, str(val))
    return cmd

# Decompose HAL file into slices.  For each slice return a copy of the
# input options.  In particular, we edit the hal2maf options to correspond
# to the slice
def computeSlices(options):
    
    

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
                        "if specifed, output MAF.", type=int,
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
                        default=0)
    h2mGrp.add_argument("--length",
                        help="length of the reference genome (or sequence"
                        " if specified) to convert.  If set to 0,"
                        " the entire thing is converted", type=int,
                        default=0)
    h2mGrp.add_argument("--rootGenome", 
                        help="name of root genome (none if empty)", 
                        default=None)
    h2mGrp.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all if empty)",
                        default=None)
    h2mGrp.add_argument("--maxRefGap", 
                        help="maximum gap length in reference", type=int,
                        default=0)
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

if __name__ == "__main__":
    sys.exit(main())
