#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Compute the background mutation rate over selected genomic regions
using bed files of mutations, and bed files of target regions.  
"""
import argparse
import os
import sys
import copy
import random
import math
from collections import defaultdict
import numpy as np
import subprocess


from hal.analysis.neutralIndel.bedMutations import BedMutations
from hal.mutations.impl.halTreeMutations import runShellCommand

# pipe stdout of given command (presumably a bedtools tool) line by
# line using a generator
def scanBedOutput(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                               stderr=sys.stderr, bufsize=-1)
    for line in process.stdout:
        yield line
    output, nothing = process.communicate()
    for line in output:
        yield line
    if process.returncode != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" %
                           (command, process.returncode))

# use above method for intersectBed.  
def scanBedIntersectionMutations(bedMutations, bedSelections):
    command = "intersectBed -a %s -b %s -wa" % (bedMutations, bedSelections)
    for line in scanBedOutput(command):
        yield line

# sum all intervals in a bed file.  assume that there are no
# overlaps! 
def computeSelectionSize(bedSelectionPath):
    size = 0
    bedSelection = open(bedSelectionPath, "r")
    for line in bedSelection:
        tokens = line.split()
        if len(tokens) == 0:
            continue
        if tokens[0][0] == "#":
            continue
        assert len(tokens) > 2
        size += int(tokens[2]) - int(tokens[1])
    return size

# given a list of mutations, and a list of regions, count the
# number of mutations in the intersection
def countMutationsInOverlap(bedMutationsPath, bedSelectionPath, events):
    eventSet = set(events)
    count = 0

    for line in scanBedIntersectionMutations(bedMutationsPath,
                                             bedSelectionPath):
        tokens = line.split()
        if len(tokens) == 0:
            continue
        if tokens[0][0] == "#":
            continue
        assert len(tokens) > 3
        tag = tokens[3]
        if tag[0] == "S" and BedMutations.substitutionBedTag in eventSet or\
           tag in eventSet:
            count += 1
    return count

# compute the rate of mutations / site given a list of mutations, and
# query regions. 
def getBackgroundRate(bedMutationsPath, bedSelectionPath, events):
    size = computeSelectionSize(bedSelectionPath)
    count = countMutationsInOverlap(bedMutationsPath, bedSelectionPath, events)
    return count, size
            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("eventsBed", type=str,
                        help="bed file containing mutation events")
    parser.add_argument("selectBed", type=str,
                        help="bed file containing regions to examine")
    parser.add_argument("--events", default=
                        " ".join(BedMutations.defaultEvents),
                        type=str, help="event tags")
    
    args = parser.parse_args()
    
    events =  args.events.split()
    count, size = getBackgroundRate(args.eventsBed, args.selectBed, events)
    print("%d / %d = %f" % (count, size, float(count) / float(size)))
    
if __name__ == "__main__":
    sys.exit(main())

