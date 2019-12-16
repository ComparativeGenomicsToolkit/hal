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

from hal.analysis.neutralIndel.bedConservation import BedConservation
from hal.analysis.neutralIndel.bedMutations import BedMutations
from hal.mutations.impl.halTreeMutations import getHalTreeMutations
from hal.mutations.impl.halTreeMutations import runShellCommand
from hal.mutations.impl.halTreeMutations import getHalRootName
from hal.mutations.impl.halTreeMutations import getHalParentName
from hal.mutations.impl.halTreeMutations import getHalChildrenNames

def getHalTreeBackground(halPath, args, rootName=None):
    root = rootName
    if root is None:
        root = getHalRootName(halPath)
    for child in getHalChildrenNames(halPath, root):
        bgFile = os.path.join(args.workDir, args.backgroundBedName % child)
        if args.ar is True:
            command = "halMaskExtract %s %s --maskFile %s --extend %d --extendPct %f" % (halPath, child, bgFile, args.arExtend, args.arExtendPct)
        else:
            command = "halStats %s --bedSequences %s > %s" % (halPath, child,
                                                              bgFile)
        print(command)
        runShellCommand(command)
        getHalTreeBackground(halPath, args, child)

        
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
    parser.add_argument("--ar", action="store_true", default=False,
                        help="Select only repeatmasked regions")
    parser.add_argument("--arExtend", type=int, default=0,
                        help="Extend selected repeats by given number of bases")
    parser.add_argument("--arExtendPct", type=float, default=0.0,
                        help="Extend selected repeated regions by given"
                        " percent")
    parser.add_argument("--root", default=None, type=str, help="root")

    args = parser.parse_args()
    args.backgroundBedName = args.backgroundBedName.replace("%%", "%")

    if args.arExtend != 0 and args.arExtendPct != 0:
        raise RuntimeError("--arExtend and --arExtendPct are exclusive")
    if args.arExtend != 0 or args.arExtendPct != 0:
        args.ar = True

    if not os.path.exists(args.workDir):
        os.makedirs(args.workDir)
        
    getHalTreeBackground(args.hal, args, args.root)
    
if __name__ == "__main__":
    sys.exit(main())
