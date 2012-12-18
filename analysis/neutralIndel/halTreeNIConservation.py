#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Make histogram of inter-event distances in BED file
"""
import argparse
import os
import sys
import copy
import subprocess

from hal.analysis.neutralIndel.bedConservation import BedConservation
from hal.analysis.neutralIndel.bedMutations import BedMutations
from hal.mutations.impl.halTreeMutations import getHalTreeMutations
from hal.mutations.impl.halTreeMutations import runShellCommand
from hal.mutations.impl.halTreeMutations import getHalRootName
from hal.mutations.impl.halTreeMutations import getHalParentName
from hal.mutations.impl.halTreeMutations import getHalChildrenNames

def checkFiles(background, mutations):
    if not os.path.isfile(background):
        raise RuntimeError("Background selection file %s not found. Make "
                           "sure halTreeNIBackground has been run "
                           "and that the paths are correctly "
                           "specified" % background)
    if not os.path.isfile(mutations):
        raise RuntimeError("Mutations file %s not found.  Make sure that "
                           "halTreeMutations has been run amd that the "
                           "paths are correctly specified" % mutations)
    
def getHalTreeConservation(halPath, args, events, rootName=None):
    root = rootName
    if root is None:
        root = getHalRootName(halPath)
    for child in getHalChildrenNames(halPath, root):
        bgFile = os.path.join(args.workDir, args.backgroundBedName % child)
        muFile = os.path.join(args.workDir, args.mutationsBedName % child)
        checkFiles(bgFile, muFile)
        outPath = os.path.join(args.workDir, args.conservedBedName % child)
        outFile = open(outPath, "w")
        bc = BedConservation()
        bc.computeBackgroundRate(muFile, bgFile, events)
        bc.identifyConservedIntervals(muFile, outFile, args.pval)
        getHalTreeConservation(halPath, args, events, child)
        
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("hal", help="input hal")
    parser.add_argument("workDir", help="working dir for all bed files")
    parser.add_argument("--backgroundBedName", type=str,
                        default="%%s_bg.bed", help="Name function for "
                        "background bed files where genome name is "
                        "specified as %%s.  Computed using halTreeNIBackground")
    parser.add_argument("--mutationsBedName", type=str,
                        default="%%s.bed", help="Name function for "
                        "background bed files where genome name is "
                        "specified as %%s.  Computed using halTreeMutations")
    parser.add_argument("--conservedBedName", type=str,
                        default="%%s_cons.bed", help="Name function for output "
                        "bed files where genome name is specifed as %%s")
    parser.add_argument("--root", default=None, type=str, help="root")
    parser.add_argument("--events",
                        default="\"%s\"" % " ".join(BedMutations.defaultEvents),
                        type=str, help="event tags (must be the same as was "
                        "used for the background!!)")
    parser.add_argument("--pval", type=float, default=0.05,
                        help="max pval of conserved segment")
    
    args = parser.parse_args()
    args.backgroundBedName = args.backgroundBedName.replace("%%", "%")
    args.mutationsBedName = args.mutationsBedName.replace("%%", "%")
    args.conservedBedName = args.conservedBedName.replace("%%", "%")
    events =  args.events.split()
        
    getHalTreeConservation(args.hal, args, events, args.root)
    
if __name__ == "__main__":
    sys.exit(main())
