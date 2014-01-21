#!/usr/bin/env python
# won't work in the general case right now
# (i.e. inputs need to be single copy & positive query strand)
import sys
import os
import argparse
from sonLib import bioio
from sonLib.bioio import getTempFile
from operator import itemgetter
from collections import defaultdict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions
class Setup(Target):
    def __init__(self, args):
        Target.__init__(self)
        self.args = args

    def run(self):
        numLines = self.numLinesInFile(self.args.bedFile)
        slices = []
        sliceNum = self.args.sliceNum
        if numLines > self.args.sliceNum:
            step = numLines/self.args.sliceNum
            count = 0
            for i in xrange(sliceNum):
                slices.append((count, count + step))
                count += step
            slices[-1] = (slices[-1][0], slices[-1][1] + numLines % self.args.sliceNum)
        else:
            sliceNum = numLines
            slices = [1]*sliceNum

        sliceOutputs = []
        for i in xrange(sliceNum):
            slice = slices[i]
            sliceOut = getTempFile(rootDir=self.getGlobalTempDir())
            sliceOutputs.append(sliceOut)
            self.addChildTarget(RunContiguousRegions(self.args, slice,
                                                     sliceOut))
        self.setFollowOnTarget(WriteToOutput(self.args, sliceOutputs))

    def numLinesInFile(self, fileName):
        n = 0
        for line in open(fileName):
            n += 1
        return n

class RunContiguousRegions(Target):
    def __init__(self, args, slice, sliceOut):
        Target.__init__(self)
        self.args = args
        self.slice = slice
        self.sliceOut = sliceOut

    def run(self):
        contiguousRegions = ContiguousRegions(self.args.alignment,
                                              self.args.srcGenome,
                                              self.args.destGenome,
                                              self.args.maxGap,
                                              self.getGlobalTempDir())
        startLineNum = self.slice[0]
        endLineNum = self.slice[1]
        outFile = open(self.sliceOut, 'w')
        for line in contiguousRegions.getContiguousLines(self.args.bedFile,
                                                         startLineNum,
                                                         endLineNum):
            outFile.write(line)

class WriteToOutput(Target):
    def __init__(self, args, sliceFiles):
        Target.__init__(self)
        self.args = args
        self.sliceFiles = sliceFiles

    def run(self):
        outFile = open(self.args.outFile, 'w')
        for filename in self.sliceFiles:
            for line in open(filename):
                outFile.write(line)

class ContiguousRegions:
    def __init__(self, alignment, srcGenome, destGenome, maxGap, tempRoot):
        self.alignment = alignment
        self.srcGenome = srcGenome
        self.destGenome = destGenome
        self.maxGap = maxGap
        self.tempRoot = tempRoot

    def liftover(self, bedLine):
        tempSrc = getTempFile("ContiguousRegions.tempSrc.bed",
                                    rootDir=self.tempRoot)
        tempDest = getTempFile("ContiguousRegions.tempDest.psl",
                                     rootDir=self.tempRoot)
        open(tempSrc, 'w').write("%s\n" % bedLine)
        cmd = "halLiftover --outPSL %s %s %s %s %s" % (self.alignment,
                                                       self.srcGenome,
                                                       tempSrc,
                                                       self.destGenome,
                                                       tempDest)
        bioio.system(cmd)
        pslLines = open(tempDest).read().split("\n")
        os.remove(tempSrc)
        os.remove(tempDest)
        pslLines = map(lambda x: x.split(), pslLines)
        # dict is to keep blocks separated by target sequence
        qBlocks = defaultdict(list)
        tBlocks = defaultdict(list)
        for pslLine in pslLines:
            if pslLine == []:
                continue
            qStrand = pslLine[8][0]
            assert(qStrand == '+')
            if len(pslLine[8]) != 1:
                assert(len(pslLine[8]) == 2)
                tStrand = pslLine[8][1]
            else:
                tStrand = '+'
            tName = pslLine[13]
            tSize = int(pslLine[14])
            blockSizes = [int(i) for i in pslLine[18].split(",") if i != '']
            qStarts = [int(i) for i in pslLine[19].split(",") if i != '']
            tStarts = [int(i) for i in pslLine[20].split(",") if i != '']
            assert(len(blockSizes) == len(qStarts) and
                   len(qStarts) == len(tStarts))
            for blockLen, qStart, tStart in zip(blockSizes, qStarts, tStarts):
                qBlocks[tName].append((qStart, qStart + blockLen))
                if tStrand == '+':
                    tBlocks[tName].append((tStart, tStart + blockLen))
                else:
                    tBlocks[tName].append((tSize - tStart - blockLen,
                                           tSize - tStart))

        # take only the blocks from the target sequence with the most mapped
        # bases
        tSeqMapped = []
        for seq, blocks in tBlocks.items():
            mappedBlocks = reduce(lambda r, v: r + (v[1] - v[0]), blocks, 0)
            tSeqMapped.append((seq, mappedBlocks))
        if len(tSeqMapped) == 0:
            # can happen if the sequence doesn't map to the target at all
            return (None, None)
        tSeqName = max(tSeqMapped, key=itemgetter(1))[0]
        qBlocks = qBlocks[tSeqName]
        tBlocks = tBlocks[tSeqName]
        
        qBlocks = self.mergeBlocks(qBlocks)
        tBlocks = self.mergeBlocks(tBlocks)
        return (qBlocks, tBlocks)

    def mergeBlocks(self, blocks):
        blocks.sort(key=itemgetter(0))
        merged = []
        prev = None
        for block in blocks:
            if prev is not None:
                # haven't thought it through yet so restrict to
                # single-copy to be safe
                assert(prev[1] <= block[0])
                if prev[1] == block[0]:
                    prev = (prev[0], block[1])
                    continue # avoid prev getting overwritten
                else:
                    merged.append(prev)
            prev = block
        
        merged.append(prev)
        return merged

    def isContiguousInTarget(self, bedLine):
        (qBlocks, tBlocks) = self.liftover(bedLine)
        if qBlocks is None or tBlocks is None:
            return False
        bedStart = int(bedLine.split()[1])
        bedEnd = int(bedLine.split()[2])

        qGaps = []
        prevqEnd = bedStart
        for block in qBlocks:
            gap = block[0] - prevqEnd
            assert(gap >= 0)
            qGaps.append(gap)
            prevqEnd = block[1]
        if bedEnd > prevqEnd:
            qGaps.append(bedEnd - prevqEnd)
        
        tGaps = []
        prevtEnd = None
        for block in tBlocks:
            if prevtEnd is not None:
                gap = block[0] - prevtEnd
                assert(gap >= 0)
                tGaps.append(gap)
            prevtEnd = block[1]

        if len([i for i in qGaps if i > self.maxGap]) > 0 or len([i for i in tGaps if i > self.maxGap]) > 0:
            return False

        return True

    def getContiguousLines(self, bedPath, startLineNum=0, endLineNum=-1):
        for lineNum, line in enumerate(open(bedPath)):
            if lineNum < startLineNum:
                continue
            elif lineNum >= endLineNum:
                break
            # can't handle bed12
            assert(len(line.split()) >= 3 and len(line.split()) < 12)
            if self.isContiguousInTarget(line):
                yield line
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="HAL alignment file")
    parser.add_argument("srcGenome", help="Reference genome")
    parser.add_argument("bedFile", help="Bed file (in ref coordinates)")
    parser.add_argument("destGenome", help="Genome to check contiguity in")
    parser.add_argument("outFile", help="Output BED file")
    parser.add_argument("--maxGap", help="maximum gap size to accept", 
                        default=10, type=int)
    parser.add_argument("--sliceNum", help="number of slices to create",
                        type=int, default=1)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)
    result = Stack(Setup(args)).startJobTree(args)
    if result:
        raise RuntimeError("Jobtree has failed jobs.")

    return 0
    
if __name__ == '__main__':
    from hal.analysis.halContiguousRegions import *
    sys.exit(main())
