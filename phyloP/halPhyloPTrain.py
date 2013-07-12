#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Attempt to automate neutral model estimation from 4d sitesfor phyloP.  
"""
import argparse
import os
import sys
import copy
import subprocess
import time
import math
import glob
from collections import defaultdict
from multiprocessing import Pool

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import runParallelShellCommands
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments
from hal.stats.halStats import getHalStats
from hal.stats.halStats import getHalTree
from hal.stats.halStats import getHalBaseComposition

    
# it seems that msa view doesn't like the second line of MAF headers
# (reads the tree as a maf block then spits an error).  so we use
# this to remove 2nd lines from generated mafs. 
def remove2ndLine(path):
    runShellCommand("sed -i -e 2d %s" % path)

def extractGeneMAFs(options):
    runShellCommand("rm -f %s" % options.outMafAllPaths)

    for bedFile in options.bedFiles:
        bedFile4d = (os.path.splitext(options.outMafPath)[0] + "_" + 
                     os.path.splitext(os.path.basename(bedFile))[0] + 
                     "4d.bed")
        if not options.no4d:
            runShellCommand("hal4dExtract %s %s %s %s" % (
                options.hal, options.refGenome, bedFile, bedFile4d))
        else:
            runShellCommand("cp %s %s" % (bedFile, bedFile4d))
            
        outMaf = (os.path.splitext(options.outMafPath)[0] + "_" +
                  os.path.splitext(os.path.basename(bedFile))[0] + ".maf")
        h2mFlags = "--ucscNames --splitBySequence --noDupes"
        if options.noAncestors is True:
            h2mFlags += " --noAncestors"
        runShellCommand("hal2mafMP.py %s %s %s "
                        "--numProc %d --refTargets %s --refGenome %s "
                        % (options.hal, outMaf, h2mFlags,options.numProc,
                           bedFile4d, options.refGenome))
        os.remove(bedFile4d)
            
    for mafFile in glob.glob(options.outMafAllPaths):
        if os.path.getsize(mafFile) < 5:
            os.remove(mafFile)
        else:
            remove2ndLine(mafFile)
            runShellCommand("msa_view -o SS --in-format MAF %s > %s" % (
                mafFile, mafFile.replace(".maf", ".SS")))
    if len(glob.glob(options.outMafAllPaths)) < 1:
        raise RuntimeError("Given BED files do not overlap alignment")

def computeMAFStats(options):
    # make one big maf
    first = True
    for mafFile in glob.glob(options.outMafAllPaths):
        if first is True:
            first = False
            runShellCommand("mv %s %s" % (mafFile, options.outMafPath))
        else:
            with open(options.outMafPath, "a") as outMaf:
                with open(mafFile, "r") as inMaf:
                    for line in inMaf:
                        l = line.lstrip()
                        if len(l) > 0 and l[0] != "#":
                            outMaf.write(line + "\n")            
            
    runShellCommand("msa_view -o SS --in-format MAF %s > %s" % (
        options.outMafPath, options.outMafSS))
    runShellCommand("rm -f %s %s" % (options.outMafAllPaths, 
                                     options.outMafPath))

# How I understand it should be run from Melissa's example.  But can't get it
# to work without msa_view to crash with a:
# msa_view(75116) malloc: *** mmap(size=18446744056529682432) failed (error code=12)
# *** error: can't allocate region
def computeAgMAFStats(options):
    halSpecies = ",".join(options.halGenomes)
    
    runShellCommand("msa_view -o SS -z --in-format SS --aggregate %s %s > %s" % (
        halSpecies, options.outMafAllPaths.replace(".maf", ".SS"), 
        options.outMafSS))
    runShellCommand("rm -f %s" % options.outMafAllPaths)
    runShellCommand("rm -f %s" % options.outMafAllPaths.replace(".maf", ".SS"))


def computeFit(options):
    runShellCommand("phyloFit --tree \"%s\" --subst-mod SSREV --sym-freqs %s --out-root %s" % (
        getHalTree(options.hal), options.outMafSS, 
        os.path.splitext(options.outMod)[0]))
    runShellCommand("rm -f %s" % options.outMafSS)

def modFreqs(options):
    baseComp = getHalBaseComposition(options.hal, options.refGenome, 1)
    runShellCommand("mv %s %s_temp" % (options.outMod, options.outMod))
    runShellCommand("modFreqs %s_temp %f %f %f %f > %s" % (options.outMod,
                                                           baseComp[0],
                                                           baseComp[1],
                                                           baseComp[2],
                                                           baseComp[3],
                                                           options.outMod))
    runShellCommand("rm -f %s_temp" % options.outMod)
        

def computeModel(options):
    runShellCommand("rm -f %s" % options.outMafAllPaths)
    extractGeneMAFs(options)
    computeAgMAFStats(options)
    computeFit(options)
    modFreqs(options)
    runShellCommand("rm -f %s" % options.outMafAllPaths)
      
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compute a neutral substitution model for use with "
        "phyloP or halPhlyoP")

    parser.add_argument("hal", help="input hal")   
    parser.add_argument("refGenome", help="Name of reference genome")
    parser.add_argument("bedDir", help="BED file or directory containing BED "
                        "files.  By "
                        "default, these files are interpreted to contain only"
                        " coordinates of coding exons, and fourfold degenerate"
                        " sites will automatically be extracted from them."
                        " To disable this behaviour and train on the entire "
                        " file, use the --no4d option.", default=None)
    parser.add_argument("outMod", help="Path to output model file")
    parser.add_argument("--no4d", help="Do not extract fourfold degenerate"
                        " positions from the input bed files.  Rather use "
                        "all bases they contain.",
                        default=False, action="store_true")
    parser.add_argument("--numProc",
                        help="Maximum number of processes for hal2maf.",
                        type=int, default=1)
    parser.add_argument("--noAncestors",
                        help="Don't write ancestral genomes in hal2maf",
                        action="store_true", default=False)

    
    args = parser.parse_args()

    if not os.path.isfile(args.hal):
        raise RuntimeError("Input hal file %s not found" % args.hal)
    if not os.path.exists(args.bedDir):
        raise RuntimeError("% not found" % args.bedDir)
    if os.path.isdir(args.bedDir):
        args.bedFiles = [os.path.join(args.bedDir, f) for f
                         in os.listdir(args.bedDir) 
                         if os.path.isfile(os.path.join(args.bedDir, f))]
    else:
        args.bedFiles = [args.bedDir]
    outTest = open(args.outMod, "w")
    if not outTest:
        raise RuntimeError("Unable to open output %s" % args.outMod)
    args.halGenomes = getHalGenomes(args.hal)
    if not args.refGenome in args.halGenomes:
        raise RuntimeError("Reference genome %s not found." % args.refGenome)
    if os.path.splitext(args.outMod)[1] != ".mod":
        raise RuntimeError("Output model must have .mod extension")

    args.outDir = os.path.dirname(args.outMod)
    args.outName = os.path.splitext(os.path.basename(args.outMod))[0]
    args.outMafName = args.outName + "_halPhyloPTrain_temp.maf"
    args.outMafPath = os.path.join(args.outDir, args.outMafName)
    args.outMafAllPaths = args.outMafPath.replace("_halPhyloPTrain_temp.maf", 
                                                  "_halPhyloPTrain_temp*.maf")
    args.outMafSS = args.outMafPath.replace("_halPhyloPTrain_temp.maf", 
                                            "_halPhyloPTrain_temp.ss")
    computeModel(args)
    
if __name__ == "__main__":
    sys.exit(main())
