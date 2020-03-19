#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Little python halStats wrapper.
"""
import argparse
import os
import sys
import copy
import subprocess
from multiprocessing import Pool



def runShellCommand(command, ascii=True):
    try:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,
                                   stderr=sys.stderr, bufsize=-1)
        output, nothing = process.communicate()
        sts = process.wait()
        if sts != 0:
            raise RuntimeError("Command: %s exited with non-zero status %i" %
                               (command, sts))
        return output.decode() if ascii else output
    except KeyboardInterrupt:
        raise RuntimeError("Aborting %s" % command)

def runParallelShellCommands(cmdList, numProc):
    if numProc == 1 or len(cmdList) == 1:
        for cmd in cmdList:
            runShellCommand(cmd)
    elif len(cmdList) > 0:
        mpPool = Pool(processes=min(numProc, len(cmdList)))
        result = mpPool.map_async(runShellCommand, cmdList)
        # specifying a timeout allows keyboard interrupts to work?!
        # http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        try:
            result.get(20000000)
        except KeyboardInterrupt:
            mpPool.terminate()
            raise RuntimeError("Keyboard interrupt")
        if not result.successful():
            raise "One or more of commands %s failed" % str(cmdList)

def getHalGenomes(halPath):
    return runShellCommand("halStats %s --genomes" % halPath).split()

def getHalNumSegments(halPath, genomeName):
    res = runShellCommand("halStats %s --numSegments %s" %
                          (halPath, genomeName)).split()
    return tuple([int(x) for x in res])

def getHalStats(halPath):
    res = runShellCommand("halStats %s" % (halPath)).split("\n")
    outList = []
    foundHeader = False
    for line in res:
        tokens = line.strip().split(",")
        if len(tokens) == 6 and tokens[0] == "GenomeName":
            foundHeader = True
        elif len(tokens) == 6 and foundHeader:
            outList.append(tuple(tokens))
    return outList

def getHalTotalStats(halPath):
    allStats = getHalStats(halPath)
    totals = (0, 0, 0, 0, 0)
    for genome in allStats:
        assert len(genome) == 6
        totals = tuple(sum(t) for t in zip(totals, genome)[1:])
    return totals

def getHalSequenceStats(halPath, genomeName):
    res = runShellCommand("halStats %s --sequenceStats %s" %
                          (halPath, genomeName)).split("\n")
    outList = []
    for line in res[1:]:
        tokens = line.strip().split(",")
        if len(tokens) == 4:
            outList.append((tokens[0], int(tokens[1]), int(tokens[2]),
                            int(tokens[3])))
    return outList

def getHalRootName(halPath):
    return runShellCommand("halStats %s --root" % halPath).strip()

def getHalParentName(halPath, genomeName):
    res = runShellCommand("halStats %s --parent %s" % (halPath, genomeName))
    return res.strip()

def getHalChildrenNames(halPath, genomeName):
    return runShellCommand("halStats %s --children %s" %
                           (halPath, genomeName)).split()

def getHalGenomeLength(halPath, genomeName):
    for genomeStats in getHalStats(halPath):
        if genomeStats[0] == genomeName:
            return int(genomeStats[2])
    return None

def getHalTree(halPath):
    return runShellCommand("halStats %s --tree" % halPath).strip()

def getHalBaseComposition(halPath, genomeName, step):
    strList = runShellCommand("halStats %s --baseComp %s,%d" % (
        halPath, genomeName, step)).split()
    return [float(x) for x in strList]

def getHalGenomeMetaData(halPath, genomeName):
    res = runShellCommand("halStats %s --genomeMetaData %s" % (halPath,
                                                               genomeName))
    if res.strip() == '':
        return dict()
    return dict([line.split("\t") for line in res.strip().split("\n")])
