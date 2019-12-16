#!/usr/bin/env python3
# won't work in the general case right now
# (i.e. needs to be positive query strand)
import sys
import os
import argparse
import itertools
from sonLib import bioio
from sonLib.bioio import getTempFile
from operator import itemgetter
from collections import defaultdict
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions
from functools import reduce

# Useful itertools recipe
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

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
            for i in range(sliceNum):
                slices.append((count, count + step))
                count += step
            slices[-1] = (slices[-1][0], slices[-1][1] + numLines % self.args.sliceNum)
        else:
            sliceNum = numLines
            slices = [[0, numLines]]

        sliceOutputs = []
        for i in range(sliceNum):
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
                                              self.getGlobalTempDir(),
                                              self.args.maxIntronDiff,
                                              self.args.deletionGaps,
                                              self.args.requiredMapFraction)
        startLineNum = self.slice[0]
        endLineNum = self.slice[1]
        outFile = open(self.sliceOut, 'w')
        for (metCriteria, line, numPreservedBases, numMapping, length) in contiguousRegions.getContiguousLines(self.args.bedFile,
                                                         startLineNum,
                                                         endLineNum):
            if self.args.printStats:
                outFile.write("%d\t%d\t%d\n" % (numPreservedBases, numMapping, length))
            elif metCriteria:
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
    def __init__(self, alignment, srcGenome, destGenome, maxGap, tempRoot,
                 maxIntronDiff, noDeletions, requiredMapFraction):
        self.alignment = alignment
        self.srcGenome = srcGenome
        self.destGenome = destGenome
        self.maxGap = maxGap
        self.tempRoot = tempRoot
        self.maxIntronDiff = maxIntronDiff
        self.noDeletions = noDeletions
        self.requiredMapFraction = requiredMapFraction

    def liftover(self, bedLine):
        """Lift a bedLine over to the target genome, parse the PSL output, and
        return a map from target sequence -> [(query block, [target
        block(s)])]

        Blocks are (start, end, strand) where start < end

        """
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
        pslLines = [x.split() for x in pslLines]
        # Get target blocks for every query block. All adjacencies
        # within a block are by definition preserved. Adjacencies
        # between target blocks (and query blocks with the commandline
        # option) are what determine if the structure is preserved.
        # dict is to keep blocks separated by target sequence & strand
        blocks = defaultdict(list)
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
                qBlock = (qStart, qStart + blockLen, qStrand)
                tBlock = (tStart, tStart + blockLen, tStrand) if tStrand == '+' else (tSize - tStart - blockLen, tSize - tStart, tStrand)
                blocks[tName].append((qBlock, tBlock))

        # Sort & merge query blocks in cases of duplication
        return self.mergeBlocks(blocks)

    def mergeBlock(self, block1, block2):
        """block1 is (qBlock, [tBlock1, tBlock2, ...]), block2 is (qBlock,
        [tBlock1, tBlock2, ...]
        """
        def takeFirst(len, block):
            if block[2] == '+':
                return (block[0], block[0] + len, block[2])
            else:
                return (block[1] - len, block[1], block[2])
        def takeLast(len, block):
            if block[2] == '+':
                return (block[1] - len, block[1], block[2])
            else:
                return (block[0], block[0] + len, block[2])
        preOverlapBlock = None
        overlapBlock = None
        postOverlapBlock = None
        qBlock1 = block1[0]
        qBlock2 = block2[0]
        assert(qBlock1[2] == '+')
        assert(qBlock2[2] == '+')
        tBlocks1 = block1[1]
        if not isinstance(tBlocks1, list):
            tBlocks1 = [tBlocks1]
        for tBlock in tBlocks1:
            assert((tBlock[1] - tBlock[0]) == (qBlock1[1] - qBlock1[0]))
        tBlocks2 = block2[1]
        for tBlock in tBlocks2:
            assert((tBlock[1] - tBlock[0]) == (qBlock2[1] - qBlock2[0]))
        if qBlock1[0] < qBlock2[1] and not qBlock1[1] <= qBlock2[0]:
            # overlapping query block
            assert(qBlock1[0] >= qBlock2[0])
            preOverlapSize = qBlock1[0] - qBlock2[0]
            postOverlapSize = abs(qBlock1[1] - qBlock2[1])
            if qBlock1[0] > qBlock2[0]:
                # block before overlap
                preOverlapqBlock = (qBlock2[0], qBlock1[0], qBlock2[2])
                preOverlaptBlocks = [takeFirst(preOverlapSize, x) for x in tBlocks2]
                preOverlapBlock = (preOverlapqBlock, preOverlaptBlocks)
            # overlapping block
            overlapSize = abs(min(qBlock1[1], qBlock2[1]) - qBlock1[0])
            if qBlock1[1] > qBlock2[1]:
                overlapqBlock = (qBlock1[0], qBlock1[1] - postOverlapSize, qBlock1[2])
                overlaptBlocks = [takeLast(overlapSize, x) for x in tBlocks2] + [takeLast(overlapSize, takeFirst(overlapSize, x)) for x in tBlocks1]
                overlapBlock = (overlapqBlock, overlaptBlocks)
            else:
                overlapqBlock = (qBlock1[0], qBlock2[1] - postOverlapSize, qBlock1[2])
                overlaptBlocks = [takeLast(overlapSize, takeFirst(preOverlapSize + overlapSize, x)) for x in tBlocks2] + tBlocks1
                overlapBlock = (overlapqBlock, overlaptBlocks)
            if qBlock1[1] > qBlock2[1]:
                # block after overlap
                postOverlapqBlock = (qBlock2[1], qBlock1[1], qBlock1[2])
                postOverlaptBlocks = [takeLast(postOverlapSize, x) for x in tBlocks1]
                postOverlapBlock = (postOverlapqBlock, postOverlaptBlocks)
            elif qBlock1[1] < qBlock2[1]:
                # block after overlap
                postOverlapqBlock = (qBlock1[1], qBlock2[1], qBlock1[2])
                postOverlaptBlocks = [takeLast(postOverlapSize, x) for x in tBlocks2]
                postOverlapBlock = (postOverlapqBlock, postOverlaptBlocks)
        else:
            preOverlapBlock = block2
            postOverlapBlock = (qBlock1, tBlocks1)
        return (preOverlapBlock, overlapBlock, postOverlapBlock)

    def insertIntoBlockList(self, blocks, blockList):
        retList = []
        if not isinstance(blocks[1], list):
            blocks = (blocks[0], [blocks[1]])
        for i, listBlock in enumerate(blockList):
            listqBlock = listBlock[0]
            qBlock = blocks[0]
            if listqBlock[1] < qBlock[0]:
                retList.append(listBlock)
                continue
            preOverlap, overlap, postOverlap = self.mergeBlock(blocks, listBlock)
            restOfList = blockList[:i] + blockList[i+1:]
            if len(restOfList) >= 1 and postOverlap is not None:
                retList = self.insertIntoBlockList(postOverlap, restOfList)
            elif postOverlap is not None:
                retList.append(postOverlap)
            if preOverlap is not None:
                retList.append(preOverlap)
            if overlap is not None:
                retList.append(overlap)
            return sorted(retList, key=itemgetter(0))
        # if we've gotten here there is no block to merge with
        retList.append((blocks[0], blocks[1]))
        return sorted(retList, key=itemgetter(0))

    def mergeBlocks(self, blockDict):
        """Take a dict of lists of (query block, target block) and turn it
        into a dict of lists of (query block, [target block(s)]),
        sorted by query block start.
        """

        ret = {}
        for seq, blockList in list(blockDict.items()):
            blockList.sort(key=itemgetter(0))
            newBlockList = []
            prev = None
            for blocks in blockList:
                if prev is not None:
                    newBlockList = self.insertIntoBlockList(blocks, newBlockList)
                else:
                    # sloppy
                    newBlockList.append((blocks[0], [blocks[1]]))
                prev = newBlockList[-1]
            ret[seq] = newBlockList
        return ret

    def isPreserved(self, blocks1, blocks2, maxGap=None, minGap=0):
        """Check if any possible adjacency between the target blocks is
           preserved. Query start for blocks1 should be less than or
           equal to query start for blocks2.
        """
        if maxGap is None:
            maxGap = self.maxGap
        for x, y in itertools.product(blocks1, blocks2):
            if x[2] == y[2]: # orientation preserved
                if x[2] == '+' and y[0] - x[1] in range(minGap, maxGap):
                    return True
                elif x[2] == '-' and x[0] - y[1] in range(minGap, maxGap):
                    return True
        return False

    def isContiguousInTarget(self, bedLine):
        elementIsPreserved = False
        blockDict = self.liftover(bedLine)
        if blockDict is None:
            return (False, 0)

        bedFields = bedLine.split()
        bedStart = int(bedFields[1])
        bedEnd = int(bedFields[2])
        bedIntrons = []
        bedLength = 0
        if len(bedFields) == 12:
            blockStarts = list(map(int, [x for x in bedFields[11].split(",") if x != ""]))
            blockSizes = list(map(int, [x for x in bedFields[10].split(",") if x != ""]))
            assert(len(blockStarts) == len(blockSizes))
            bedBlocks = [(bedStart + start, bedStart + start + size)
                         for start, size in zip(blockStarts, blockSizes)]
            prevEnd = None
            for block in bedBlocks:
                if prevEnd is not None:
                    gap = block[0] - prevEnd
                    assert(gap >= 0)
                    bedIntrons.append((prevEnd, block[0]))
                prevEnd = block[1]
            bedLength = sum(blockSizes)
        else:
            bedLength = bedEnd - bedStart

        numPreservedAdjacencies = 0

        # take only the blocks from the target sequence with the most mapped
        # bases
        totalMappedAdjacencies = 0
        tSeqMapped = {}
        for seq, value in list(blockDict.items()):
            qBlocks = list(map(itemgetter(0), value))
            mappedBases = reduce(lambda r, v: r + (v[1] - v[0]), qBlocks, 0)
            totalMappedAdjacencies += mappedBases - 1
            # Adjacencies within blocks are always preserved.
            numPreservedAdjacencies += mappedBases - len(qBlocks)
            mappedFraction = float(mappedBases)/bedLength
            tSeqMapped[seq] = mappedFraction
        if len(tSeqMapped) == 0:
            # can happen if the sequence doesn't map to the target at all
            return (False, 0, 0, bedLength)

        for seq, blocks in list(blockDict.items()):
            # FIXME: Need to account for introns again
            # And qGaps if option is given
            preservedForSeq = True
            if tSeqMapped[seq] < self.requiredMapFraction:
                preservedForSeq = False
            for (qBlock1, tBlocks1), (qBlock2, tBlocks2) in pairwise(blocks):
                maxGap = self.maxGap
                minGap = 0
                for intron in bedIntrons:
                    assert(qBlock2[0] >= qBlock1[1])
                    qGap = qBlock2[0] - qBlock1[1]
                    if qBlock2[0] >= intron[1] and qBlock1[1] <= intron[0]:
                        # query gap is from this intron
                        maxGap = qGap + self.maxIntronDiff
                        minGap = qGap - self.maxIntronDiff
                if self.isPreserved(tBlocks1, tBlocks2, minGap=minGap,
                                    maxGap=maxGap):
                    numPreservedAdjacencies += 1
                else:
                    preservedForSeq = False
            if preservedForSeq:
                elementIsPreserved = True

        return (elementIsPreserved, numPreservedAdjacencies, totalMappedAdjacencies, bedLength)

    def getContiguousLines(self, bedPath, startLineNum=0, endLineNum=-1):
        for lineNum, line in enumerate(open(bedPath)):
            if lineNum < startLineNum:
                continue
            elif lineNum >= endLineNum:
                break

            (metCriteria, numAdjacencies, numMapped, length) = self.isContiguousInTarget(line)
            yield (metCriteria, line, numAdjacencies, numMapped, length)
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="HAL alignment file")
    parser.add_argument("srcGenome", help="Reference genome")
    parser.add_argument("bedFile", help="Bed file (in ref coordinates)")
    parser.add_argument("destGenome", help="Genome to check contiguity in")
    parser.add_argument("outFile", help="Output BED file")
    parser.add_argument("--maxGap", help="maximum gap size to accept", 
                        default=100, type=int)
    parser.add_argument("--deletionGaps", help="care about deletion gaps",
                        default=False, action='store_true')
    parser.add_argument("--sliceNum", help="number of slices to create",
                        type=int, default=1)
    parser.add_argument("--maxIntronDiff", help="Maximum number of bases "
                        "that intron gaps are allowed to change by", type=int,
                        default=10000)
    parser.add_argument("--requiredMapFraction", help="Fraction of bases in "
                        "the query that need to map to the target to be "
                        "accepted", type=float, default=0.0)
    parser.add_argument("--printStats", help="instead of printing the "
                        "passing BED lines, print statistics",
                        action='store_true', default=False)
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
