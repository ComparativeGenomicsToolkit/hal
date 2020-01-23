#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Simulate snake track browser queries and time them.  Ideally done with hgTracks
but it's easier to automate this way, as far as I can tell. 
"""
import argparse
import os
import sys
import copy
import subprocess
import time
import math
import random
from datetime import datetime
from collections import OrderedDict

from hal.stats.halStats import runShellCommand
from hal.stats.halStats import getHalGenomes
from hal.stats.halStats import getHalNumSegments
from hal.stats.halStats import getHalStats


# Wrapper for blockVizTime
def getBlockVizCmd(options, tgtGenome):
    cmd = "./blockVizTime %s %s %s %s %d %d" % (options.lod, tgtGenome,
                                              options.refGenome, options.refSequence,
                                              options.refFirst, options.refLast)
    if options.refLength >= options.maxSnp:
        cmd += " 1"
    else:
        cmd += " 0"
    if options.doDupes is True:
        cmd += " 1"
    else:
        cmd += " 0"
    if options.udc is not None:
        cmd += " %s" % options.udc

    return cmd

def timeCmd(cmd):
    tstart = datetime.now()
    runShellCommand(cmd)
    tend = datetime.now()
    tdelta = tend - tstart
    tsecs = tdelta.seconds
    tsecs += tdelta.microseconds / 1000000.0
    return tsecs

def simulateLoad(options):
    cmds = [getBlockVizCmd(options, tgtGenome) for tgtGenome in options.tgtGenomes]
    elapsedTime = 0.
    for cmd in cmds:
        lastExcep = None
        for trial in range(options.retry):
            if options.udc is not None and options.zapUdc is True:
                runShellCommand("rm -rf %s" % os.path.join(options.udc, "*"))
            t = -1
            try:
                t = timeCmd(cmd)
                lastExcep = None
                break
            except Exception as e:
                lastExcep = e
                time.sleep(2)
        if lastExcep is None:
            elapsedTime += t
        else:
            raise lastExcep
    return elapsedTime

def uniformRuns(options):
    trials = []
    for i in range(options.reps):
        size = random.uniform(1, options.refLength)
        start = random.uniform(0, options.refLength - size)
        trials.append((int(start), int(start + size)))
    return trials

def runSim(options):
    trials = uniformRuns(options)
    elapsedTimes = []
    for trial in trials:
        options.refFirst = trial[0]
        options.refLast = trial[1]
        elapsedTimes.append((options.refLast - options.refFirst, simulateLoad(options)))
    return elapsedTimes

def binTimes(options, elapsedTimes):
    binnedTimes = dict()
    for (size, time) in elapsedTimes:
        theBin = int(size / options.binSize) * options.binSize
        if theBin not in binnedTimes:
            binnedTimes[theBin] = (1, size, time, time, time)
        else:
            curVal = binnedTimes[theBin]
            binnedTimes[theBin] = (curVal[0] + 1,
                                   curVal[1] + size,
                                   curVal[2] + time,
                                   min(curVal[3], time),
                                   max(curVal[4], time))
    return binnedTimes

def printTable(options, binnedTimes):
    print("Bin, nElem, totTime, minTime, maxTime, avgTime")
    for key in sorted(binnedTimes.keys()):
        theBin = key
        nElem, totSize, totTime, minTime, maxTime = binnedTimes[theBin]
        avgTime = float(totTime) / nElem
        print("%d, %d, %f, %f, %f, %f" % (
            theBin, nElem, totTime, minTime, maxTime, avgTime))
        
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Time some simulated browser queries. ")

    parser.add_argument("lod", help="input lod or hal path")
    parser.add_argument("refGenome", help="Name of reference genome")
    parser.add_argument("refSequence", help="Name of reference sequence")
    parser.add_argument("refLength", help="Length of reference sequence", type=int)
    parser.add_argument("tgtGenomes", help="Comma-separated list of target genomes")

    parser.add_argument("--reps", help="Number of queries to perform",
                        type=int, default=100)

    parser.add_argument("--udc", help="UDC path", default=None)

    parser.add_argument("--zapUdc", help="If a UDC path is specified with --udc, then"
                        " erase it before each individual query to make sure it doesnt"
                        " get used.  this is for comparison only", action="store_true",
                        default=False)
    
    parser.add_argument("--doDupes", help="Do duplications", action="store_true",
                        default=False)

    parser.add_argument("--maxSnp", help="Max query size to get bases",
                        type=int, default=50000)

    parser.add_argument("--binSize", help="Bin size for output table",
                        type=int, default=1000)

    parser.add_argument("--seed", help="Random seed", type=int, default=None)

    parser.add_argument("--retry", help="Number of retries if a query fails (which"
                        " happens more than you think", type=int, default=10)

    args = parser.parse_args()
    args.tgtGenomes = args.tgtGenomes.split(",")

    if args.udc is not None and not os.path.exists(args.udc):
        os.makedirs(args.udc)

    if args.seed is not None:
        random.seed(args.seed)

    times = runSim(args)
    binnedTimes = binTimes(args, times)
    printTable(args, binnedTimes)
    
if __name__ == "__main__":
    sys.exit(main())
