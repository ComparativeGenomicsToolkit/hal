#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Make histogram of inter-event distances in BED file
"""
import argparse
import os
import sys
import copy
import subprocess

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import getHalRootName
from hal.stats.halStats import getHalParentName
from hal.stats.halStats import getHalChildrenNames

                        
def getHalBranchMutations(halPath, genomeName, args):
    command = "halBranchMutations %s %s --maxGap %s" % (halPath, genomeName,
                                                        args.maxGap)
    
    refBedFile = os.path.join(args.outDir,  "%s.bed" % genomeName)
    dest = refBedFile
    if not args.noSort:
        dest = "stdout"
        
    command += " --refFile %s" % dest
    command += " --delBreakFile %s" % dest
    if args.doSnps:
        command += " --snpFile %s" % dest
    if args.doParentDeletions:
        command += " --parentFile %s" % os.path.join(args.outDir, 
                                                     "%s_pd.bed" % genomeName)

    if not args.noSort:
        command += " | sortBed > %s" % refBedFile
    print(command)
    runShellCommand(command)

def getHalTreeMutations(halPath, args, rootName=None):
    root = rootName
    if root is None:
        root = getHalRootName(halPath)
    for child in getHalChildrenNames(halPath, root):
        getHalBranchMutations(halPath, child, args)
        getHalTreeMutations(halPath, args, child)

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("hal", help="input hal")
    parser.add_argument("outDir", help="output dir")
    parser.add_argument("--bedName", type=str,
                        default="%%s.bed", help="Name function for output "
                        "bed files where sequence name is specifed as %%s")

    parser.add_argument("--root", default=None, type=str, help="root")
    parser.add_argument("--doSnps",action="store_true", default=False)
    parser.add_argument("--doParentDeletions",action="store_true",
                        default=False)
    parser.add_argument("--maxGap", default=10, type=int, help="gap threshold")
    parser.add_argument("--noSort", action="store_true", default=False)
    args = parser.parse_args()

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    if not args.noSort:
        try:
            runShellCommand("echo \"x\t0\t1\" | sortBed 2> /dev/null")
        except Exception:
            print(("Warning: output BED files not sorted because sortBed" + 
               " (BedTools) not found"))
            args.noSort = True
        
    getHalTreeMutations(args.hal, args, args.root)
    
if __name__ == "__main__":
    sys.exit(main())
