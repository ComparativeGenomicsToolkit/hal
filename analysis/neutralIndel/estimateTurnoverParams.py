#!/usr/bin/env python

#Copyright (C) 2013 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Estimate constraint turnover parameters from output of halTreeNITurnover.
This file has form like
droGri2: cons 13898234  ucons 56436697  gain 11822354 (0.173198) loss 30226884 (0.685027) bl 0.183954
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
import tempfile

from hal.mutations.impl.halTreeMutations import getHalRootName
from hal.mutations.impl.halTreeMutations import getHalParentName
from hal.mutations.impl.halTreeMutations import getHalChildrenNames
from hal.anaylsis.constraintTurnover.turnoverModel import gradDescent

# go through the halTreENITurnover output, mapping observed values to genome
# in convention used by turnover model, ie
# ([pi0,pi1], [[P00, P01],[P10, P11]], t)
def readTurnoverFile(turnoverPath):
    result = dict()
    toFile = open(turnoverPath, "r")
    for line in toFile:
        toks = line.split()
        genome = toks[0].strip(":")
        cons = toks[2]
        ucons = toks[4]
        pi0 = float(ucons) / float(cons + ucons)
        pi1 = float(cons) / float(cons + ucons)
        pg = toks[7].strip("()")
        pl = toks[10].strip("()")
        t = toks[12]
        result[genome] = ([pi0, pi1], [ [1.0 - pg, pg], [pl, 1.0 - pl] ], t)
    return result

# concatenate information for all genomes below and including the given
# root (if it exists) into a list.  The observations were read from the
# output file using readTurnoverFile()
def getValuesBelowRoot(halPath, rootName, obseverations):
    nextQueue = collections.deque()
    nextQueue.append(rootName)
    output = []
    while len(nextQueue) > 0:
        next = nextQueue.popleft()
        if next in observations:
            output.append(observations[next])
        for child in getHalChildrenNames(halPath, next):
            nextQueue.append(child)
    return output

# estimate parameters using turnoverModels gradient descent function
# with random starting points between 0 and 1
# returns (rLoss, rGain, dSq)
def estimateParamsFromList(obsVals, maxIt, step, retries):
    bestLr, bestGr, bestDiff = (0, 0, 1000000)
        
    for retry in range(retries):
        lrStart = random.uniform(0.0001, 1.0)
        grStart = random.uniform(0.0001, 1.0)
        (lrEst, grEst, diff) = gradDescent(lrStart, grStart, obsVals,
                                           maxIt, step)
    if diff < bestDiff:
        bestLr, bestGr, bestDiff = (lrEst, grEst, diff)

    return [bestLr, bestGr, bestDiff]

# estimate the parameters for the root. if allInternals is true, then
# repeat for all internal nodes below the root. 
def halTreeTurnoverParams(halPath, obsPath, rootName, allInternals,
                          maxIt, step, retries):
    observations = readTurnoverFile(obsPath)
    nextQueue = collections.deque()
    nextQueue.append(rootName)
    output = []
    while len(nextQueue) > 0:
        next = nextQueue.popleft()
        if next == rootName or (allInternals is True and
                                len(getHalChildrenNames(next)) == 0):
            obsVals = getValuesBelowRoot(halPath, next, observations)
            result = estimateParamsFromList(obsVals, maxIt, step, retries)
            print "%s: lr=%f gr=%f dsq=%f" % [next] + result
        for child in getHalChildrenNames(halPath, next):
            nextQueue.append(child)
    return output
                                                              
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("halFile", sty=str,
                        help="Path of hal file")
    parser.add_argument("NITurnoverFile", type=str,
                        help="Output of halTreeNITurnover.py")
    parser.add_argument("--maxIt", type=int, default=1000,
                        help="number of iterations for gradient descent")
    parser.add_argument("--step", type=float, default=0.001,
                        help="gradient descent step")
    parser.add_argument("--retries", type=int, default=5,
                        help="number of gradient descents to run")
    parser.add_argument("--root", type=str, default=None,
                        help="root of alignment to consder")
    parser.add_argument("--allInteranls", type=store_true, default=False,
                        help="estimate params for all subtrees independently,"
                        " in addition to the root")
                       
    args = parser.parse_args()

    if args.root is None:
        args.root = getHalRootName(args.halFile)

    assert (args.maxIt > 0 and args.step > 0 and args.noise >= 0 and
            args.retries > 1)

    halTreeTurnoverParams(args.halFile, args.NITurnoverFile,
                          args.root, args.allInternals, args.maxIt,
                          args.step, args.retries)
    
  
if __name__ == "__main__":
    sys.exit(main())

