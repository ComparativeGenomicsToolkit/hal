#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Write conserved regions based on neutral-indel-like test
"""
import argparse
import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np

from hal.analysis.neutralIndel.bedMutations import BedMutations
from hal.analysis.neutralIndel.backgroundRate import getBackgroundRate
from hal.mutations.impl.halTreeMutations import runShellCommand

class BedConservation(object):

    def __init__(self):
        pass

    # use the BedHistogram code to load in a bed file and compute a
    # length distribution
    def computeBackgroundRate(self, mutationsBed, backgroundBed, events):
        self.count, self.size = getBackgroundRate(mutationsBed,
                                                  backgroundBed, events)
        
        self.rate = float(self.count) / float(self.size)
        self.events = events

    def bfProb(self, distance):
        assert self.count > 0
        p = self.rate * math.pow(1. - self.rate, distance)
        #bonferroni
        p *= self.count
        return p
        
    # after the background rate has been computed, we run a test on
    # a bed and write out conserved intervals
    def identifyConservedIntervals(self, bedPath, outStream, maxPVal=0.05):
        assert self.rate is not None
        self.writtenCount = 0
        self.writtenBases = 0
        bm = BedMutations()
        for line in bm.scan(bedPath, self.events):
            d = bm.distance()
            if d is not None:
                pval = self.bfProb(d)
                if pval <= maxPVal:
                    outStream.write("%s\t%s\t%s\t%f\t%s\t%s\n" % (
                        bm.sequence, bm.prevRange[1], bm.range[0],
                        pval, bm.ancGenome, bm.genome))
                    self.writtenBases += d
                    self.writtenCount += 1                

    def minDistance(self, maxPVal):
        i = 1
        while i < self.size:
            if self.bfProb(i) <= maxPVal:
                return i
            i += 1
        return None
            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("mutationsBed", help="input bed")
    parser.add_argument("backgroundBed", help="background regions")
    parser.add_argument("--outBed", type=str, default=None,help="output bed")
    parser.add_argument("--events",
                        default="\"%s\"" % " ".join(BedMutations.defaultEvents),
                        type=str, help="event tags")
    parser.add_argument("--pval", type=float, default=0.05,
                        help="max pval of conserved segment")
    
    args = parser.parse_args()
    outStream = sys.stdout
    if args.outBed is not None:
        outStream = open(args.outBed, "w")
    events =  args.events.split()
    
    bc = BedConservation()
    bc.computeBackgroundRate(args.mutationsBed, args.backgroundBed,
                             events)
    bc.identifyConservedIntervals(args.mutationsBed, outStream, args.pval)

    sys.stderr.write("%d segments with %d bases (%f pct of genome) found."
                     " bgrate= %f minDist=%d\n" % (
                         bc.writtenCount,
                         bc.writtenBases,
                         float(bc.writtenBases) / bc.size,
                         bc.rate,
                         bc.minDistance(args.pval)))
    
    if not outStream is sys.stdout:
        outStream.close()
    
if __name__ == "__main__":
    sys.exit(main())
