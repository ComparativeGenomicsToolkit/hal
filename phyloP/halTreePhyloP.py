#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Compute phyloP scores for each genome in an alignment.  PhyloP scores
are obtained by lifting down from the parent whenever possible to avoid 
recomputing them.    
"""
import argparse
import os
import sys
import copy
import subprocess
import time
import string
import random
import math
import glob
from collections import defaultdict
from multiprocessing import Pool

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import runParallelShellCommands
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalRootName
from hal.stats.halStats import getHalStats
from hal.stats.halStats import getHalChildrenNames
from hal.stats.halStats import getHalParentName


def outFileName(args, genome, extension, token, temp):
    fileName = "%s_%s.%s" % (genome, token, extension)
    if temp is True:
        fileName = "%s_%s" % (args.tempID, fileName)
    return os.path.join(args.outWigDir, fileName)

def computeTreePhyloP(args):
    visitQueue = [args.root]
    bigwigCmds = []
    while len(visitQueue) > 0:
        genome = visitQueue.pop()
        bedFlags = ""
        # Generate a bed file of all regions of 
        # genome that dont align to parent
        bedInsertsFile = outFileName(args, genome, "bed", "inserts", True)
        if genome != args.root:
            runShellCommand(
            "halAlignedExtract %s %s --alignedFile %s --complement" % (
                args.hal, genome, bedInsertsFile))
            bedFlags = "--refBed %s" % bedInsertsFile

        # Run halPhyloP on the inserts
        wigFile = outFileName(args, genome, "wig", "phyloP", False)
        cmd = "halPhyloPMP.py %s %s %s %s --numProc %d %s" % (
            args.hal, genome, args.mod, bedFlags, args.numProc, wigFile)
        if args.subtree is not None:
            cmd += " --subtree %s" % args.subtree
        if args.prec is not None:
            cmd += " --prec %d" % args.prec

        runShellCommand(cmd)
    
        runShellCommand("rm -f %s" % bedInsertsFile)

        # Lift down from the parent, appending to the wig file computed above
        if genome != args.root:
            parent = getHalParentName(args.hal, genome)
            parentWig = outFileName(args, parent, "wig", "phyloP", False)
            if os.path.isfile(parentWig):
                runShellCommand("halWiggleLiftover %s %s %s %s %s --append" % (
                    args.hal, parent, parentWig, genome, wigFile))

        # Convert to bigwig if desired and delete wig file
        if args.bigWig is True and os.path.isfile(wigFile):
            sizesFile = outFileName(args, genome, "sizes", "chr", True)
            bwFile = outFileName(args, genome, "bw", "phyloP", False)
            bwCmd = "halStats %s --chromSizes %s > %s && " % (args.hal, genome,
                                                              sizesFile)
            bwCmd += "wigToBigWig %s %s %s && " % (wigFile, sizesFile, bwFile)
            bwCmd += "rm -f %s &&" % wigFile
            bwCmd += "rm -f %s" % sizesFile
            bigwigCmds.append(bwCmd)

        # Recurse on children.
        children = getHalChildrenNames(args.hal, genome)
        for child in children:
            visitQueue.append(child)

    #parallel bigwig conversion
    runParallelShellCommands(bigwigCmds, args.numProc)

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute PhyloP scores (in wig format) for each genome in"
        " an alignment.  Scores are computed once per column with "
        "halPhyloPMP.py and iteratively lifted down the tree using "
        "halWiggleLiftover (starting at the root).")

    parser.add_argument("hal", help="input hal")
    parser.add_argument("mod", help="model file for PhyloP.  Can be created "
                        "with halPhyloPTrain.py")
    parser.add_argument("outWigDir", help="directory where output wig files"
                        " will be written")
    parser.add_argument("--root", help="Name of root.  If not specified the"
                        " HAL root will be used", default=None)
    parser.add_argument("--numProc",
                        help="Maximum number of processes.",
                        type=int, default=1)
    parser.add_argument("--bigWig",
                        help="Run wigToBigWig on each generated wiggle",
                        action="store_true", default=False)
    parser.add_argument("--subtree",
                        help="Run clade-specific acceleration/conservation on subtree below this node",
                        default=None)
    parser.add_argument("--prec",
                        help="Number of decimal places in wig output", type=int,
                        default=None)

    # need phyloP options here:
    
    args = parser.parse_args()

    if not os.path.isfile(args.hal):
        raise RuntimeError("Input hal file %s not found" % args.hal)
    if not os.path.isfile(args.mod):
        raise RuntimeError("Input mod file %s not found" % args.mod)
    if not os.path.isdir(args.outWigDir):
        os.makedirs(args.outWigDir)
    if not os.path.isdir(args.outWigDir):
        raise RuntimeError("%s not found" % args.outWigDir)

    args.halGenomes = getHalGenomes(args.hal)
    if args.root is None:
        args.root = getHalRootName(args.hal)

    if not args.root in args.halGenomes:
        raise RuntimeError("Root genome %s not found." % args.root)

    if args.subtree is not None and args.root not in args.halGenomes:
        raise RuntimeError("Subtree root %s not found." % args.subtree)

    # make a little id tag for temporary maf slices
    S = string.ascii_uppercase + string.digits
    args.tempID = 'halTreePhyloP' + ''.join(random.choice(S) for x in range(5))

    computeTreePhyloP(args)
    
if __name__ == "__main__":
    sys.exit(main())
