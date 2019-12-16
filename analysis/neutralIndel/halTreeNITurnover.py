#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Compute constraint turnover stats over entire tree
"""
import argparse
import os
import sys
import copy
import subprocess

from hal.analysis.neutralIndel.turnoverRate import compareConservationOverBranch
from hal.analysis.neutralIndel.turnoverRate import getParentGenomeName
from hal.analysis.neutralIndel.turnoverRate import getBranchLength

from hal.mutations.impl.halTreeMutations import runShellCommand
from hal.mutations.impl.halTreeMutations import getHalRootName
from hal.mutations.impl.halTreeMutations import getHalParentName
from hal.mutations.impl.halTreeMutations import getHalChildrenNames

def checkFile(conservations):
    if not os.path.isfile(conservations):
        raise RuntimeError("Conserved intervals file %s not found. Make "
                           "sure halTreeNIConservation has been run "
                           "and that the paths are correctly "
                           "specified" % conservations)
    
def getHalTreeTurnover(halPath, args, rootName=None):
    root = rootName
    if root is None:
        root = getHalRootName(halPath)
    for child in getHalChildrenNames(halPath, root):
        if root != getHalRootName(halPath):
        
            consFile = os.path.join(args.workDir,
                                    args.conservedBedName % child)
            checkFile(consFile)
            pconsFile = os.path.join(args.workDir,
                                     args.conservedBedName % root)
            checkFile(pconsFile)

            outMappedAlignedBed = os.path.join(args.workDir,
                                               child + "_pa.bed")
            outParentSlicedBed = os.path.join(args.workDir,
                                              child + "_pslice.bed")
            outMappedGenomeBed = os.path.join(args.workDir,
                                              child + "_pm.bed")
            outConservationBed = os.path.join(args.workDir,
                                              child + "_int.bed")
            outAlignedBed = os.path.join(args.workDir, child + "_al.bed")
            outGainBed = os.path.join(args.workDir, child + "_gain.bed")
            outLossBed = os.path.join(args.workDir, child + "_loss.bed")

            (conLen, gainLen,
             lossLen, unconLen) = compareConservationOverBranch(
                halPath, child, consFile, pconsFile,
                outMappedAlignedBed, outParentSlicedBed,
                outMappedGenomeBed, outConservationBed, outAlignedBed,
                outGainBed, outLossBed)

            gainRate = 0
            if conLen + lossLen > 0:
                gainRate = float(gainLen) / (unconLen + gainLen)
            lossRate = 0
            if unconLen + gainLen > 0:
                lossRate = float(lossLen) / (conLen + lossLen)

            branchLength = getBranchLength(halPath, child)
                
            print("%s: cons %d  ucons %d  gain %d (%f) loss %d (%f) bl %f" % (
                child,                                                
                conLen,
                unconLen,
                gainLen,
                gainRate,
                lossLen,
                lossRate,
                branchLength))
        
        getHalTreeTurnover(halPath, args, child)
        
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("hal", help="input hal")
    parser.add_argument("workDir", help="working dir for all bed files")
    parser.add_argument("--conservedBedName", type=str,
                        default="%%s_cons.bed", help="Name function for output "
                        "bed files where genome name is specifed as %%s")
    parser.add_argument("--root", default=None, type=str, help="root")
    
    args = parser.parse_args()
    args.conservedBedName = args.conservedBedName.replace("%%", "%")
        
    getHalTreeTurnover(args.hal, args, args.root)
    
if __name__ == "__main__":
    sys.exit(main())
