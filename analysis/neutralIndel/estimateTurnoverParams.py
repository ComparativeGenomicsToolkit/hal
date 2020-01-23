#!/usr/bin/env python3

#Copyright (C) 2013 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

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
from collections import deque
import numpy as np
import subprocess
import tempfile

from hal.mutations.impl.halTreeMutations import getHalRootName
from hal.analysis.neutralIndel.turnoverRate import getParentGenomeName
from hal.analysis.neutralIndel.turnoverRate import getBranchLength
from hal.mutations.impl.halTreeMutations import getHalChildrenNames
from hal.analysis.constraintTurnover.turnoverModel import gradDescent
from hal.analysis.constraintTurnover.turnoverModel import computePMatrix
from hal.analysis.constraintTurnover.turnoverModel import computeStationaryDist

# go through the halTreENITurnover output, mapping observed values to genome
# in convention used by turnover model, ie
# ([pi0,pi1], [[P00, P01],[P10, P11]], t)
def readTurnoverFile(halPath, turnoverPath):
    result = dict()
    toFile = open(turnoverPath, "r")
    for line in toFile:
        toks = line.split()
        genome = toks[0].strip(":")
        cons = float(toks[2])
        ucons = float(toks[4])
        gain = float(toks[6])
        loss = float(toks[9])
        totalAligned = cons + ucons + gain + loss
                     
        if totalAligned <= 0 or cons < 0 or ucons < 0:
            sys.stderr.write("Warning, skipping %s\n" % genome)
        else:
            pi0 = (ucons + loss) / totalAligned
            pi1 = (cons + gain) / totalAligned
            pg = gain / (ucons + gain)
            pl = loss / (cons + loss)
            t = float(toks[12])
            #
            # Incorporate parent branch since it affects turnover
            #
            if genome != getHalRootName(halPath):
                parName = getParentGenomeName(halPath, genome)
                parBranch = getBranchLength(halPath, parName)
                t += float(parBranch)
            assert pi0 >= 0 and pi1 >=0
            assert pg >= 0 and pl >=0
            assert t >= 0
            result[genome] = ([pi0, pi1], [ [1.0 - pg, pg], [pl, 1.0 - pl] ], t)
    return result

# concatenate information for all genomes below and not including the given
# root (if it exists) into a list.  The observations were read from the
# output file using readTurnoverFile()
def getValuesBelowRoot(halPath, rootName, observations):
    nextQueue = deque()
    nextQueue.append(rootName)
    output = []
    while len(nextQueue) > 0:
        next = nextQueue.popleft()
        if next != rootName:
            if next in observations:
                output.append(observations[next])
            else:
                sys.stderr.write("Warning, no observation for %s\n" % next)
        for child in getHalChildrenNames(halPath, next):
            nextQueue.append(child)
    return output

# estimate parameters using turnoverModels gradient descent function
# with random starting points between 0 and 1
# returns (rLoss, rGain, dSq)
def estimateParamsFromList(obsVals, maxIt, step, retries):
    assert len(obsVals) > 0
    bestLr, bestGr, bestDiff = (0, 0, 1000000)
        
    for retry in range(retries):
        if retry == 0:
            lrStart = step
            grStart = step
        else:
            lrStart = random.uniform(0., step * maxIt)
            grStart = random.uniform(0., step * maxIt)
        (lrEst, grEst, diff) = gradDescent(lrStart, grStart, obsVals,
                                           maxIt, step)
        if diff < bestDiff:
            bestLr, bestGr, bestDiff = (lrEst, grEst, diff)

    return [bestLr, bestGr, bestDiff]

# output some stats on how the model fits the data for various nodes on the
# tree.  will need to start thinking about properly putting this in a
# spreadsheet or something
def printComparison(halPath, obsVals, observations, result):
    lossRate = result[0]
    gainRate = result[1]
    obsScope = set([str(x) for x in obsVals])
    if len(list(observations.items())) > 0:
        print("Genome, t, piObs0, piObs1, piEst0, piEst1, PLossObs, PGainObs, PLossEst, PGainEst, AvgDiff")
    for (name, obs) in list(observations.items()):
        if str(obs) in obsScope:            
            t = obs[2]
            pi = computeStationaryDist(lossRate, gainRate, t)
            P = computePMatrix(lossRate, gainRate, t)
            print("  %s, %f, %.2f, %.2f, %.2f, %.2f, %.3f, %.3f, %.3f, %.3f, %.3f" % (
                name, t, obs[0][0], obs[0][1], pi[0], pi[1], 
                obs[1][0][1], obs[1][1][0], P[0][1], P[1][0],
                0.25 * (math.fabs(obs[0][0] - pi[0]) +
                math.fabs(obs[0][1] - pi[1]) +
                math.fabs(obs[1][0][1] - P[0][1]) +
                math.fabs(obs[1][1][0] - P[1][0]))
                ))
    
# estimate the parameters for the root. if allInternals is true, then
# repeat for all internal nodes below the root. 
def halTreeTurnoverParams(halPath, obsPath, rootName, allInternals,
                          maxIt, step, retries):
    observations = readTurnoverFile(halPath, obsPath)
    nextQueue = deque()
    nextQueue.append(rootName)
    output = []
    while len(nextQueue) > 0:
        next = nextQueue.popleft()
        if next == rootName or (allInternals is True and
                                len(getHalChildrenNames(halPath, next)) > 0):

            obsVals = getValuesBelowRoot(halPath, next, observations)
            result = estimateParamsFromList(obsVals, maxIt, step, retries)
            print("%s: lr=%f gr=%f dsq=%f" % (next,result[0], result[1],
                                              result[2]))
            printComparison(halPath, obsVals, observations, result)
        for child in getHalChildrenNames(halPath, next):
            nextQueue.append(child)
    return output
                                                              
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("halFile", type=str,
                        help="Path of hal file")
    parser.add_argument("NITurnoverFile", type=str,
                        help="Output of halTreeNITurnover.py")
    parser.add_argument("--maxIt", type=int, default=100000,
                        help="number of iterations for gradient descent")
    parser.add_argument("--step", type=float, default=0.0001,
                        help="gradient descent step")
    parser.add_argument("--retries", type=int, default=5,
                        help="number of gradient descents to run")
    parser.add_argument("--root", type=str, default=None,
                        help="root of alignment to consder")
    parser.add_argument("--allInternals", action="store_true", default=False,
                        help="estimate params for all subtrees independently,"
                        " in addition to the root")
                       
    args = parser.parse_args()

    if args.root is None:
        args.root = getHalRootName(args.halFile)

    assert (args.maxIt > 0 and args.step > 0 and args.retries > 1)

    halTreeTurnoverParams(args.halFile, args.NITurnoverFile,
                          args.root, args.allInternals, args.maxIt,
                          args.step, args.retries)
    
  
if __name__ == "__main__":
    sys.exit(main())

