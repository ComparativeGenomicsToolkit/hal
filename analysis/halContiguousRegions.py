#!/usr/bin/env python
# won't work in the general case right now
# (i.e. needs to be positive query strand)
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
            slices = [[0, numLines]]

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
                                              self.getGlobalTempDir(),
                                              self.args.maxIntronDiff,
                                              self.args.noDeletions,
                                              self.args.requiredMapFraction)
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
    def __init__(self, alignment, srcGenome, destGenome, maxGap, tempRoot, maxIntronDiff, noDeletions, requiredMapFraction):
        self.alignment = alignment
        self.srcGenome = srcGenome
        self.destGenome = destGenome
        self.maxGap = maxGap
        self.tempRoot = tempRoot
        self.maxIntronDiff = maxIntronDiff
        self.noDeletions = noDeletions
        self.requiredMapFraction = requiredMapFraction

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
        f = open('log', 'a')
        pslLines = map(lambda x: x.split(), pslLines)
        tStrands = dict()
        # dict is to keep blocks separated by target sequence
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
            if tName in tStrands and tStrands[tName] != tStrand:
                # does not preserve orientation
                if tName in blocks: # there could be a duplication as well (3 lines)
                    del blocks[tName]
                continue
            tStrands[tName] = tStrand
            for blockLen, qStart, tStart in zip(blockSizes, qStarts, tStarts):
                qBlock = (qStart, qStart + blockLen)
                tBlock = (tStart, tStart + blockLen) if tStrand == '+' else (tSize - tStart - blockLen, tSize - tStart)
                blocks[tName].append((qBlock, tBlock))

        # Filter out seqs that don't preserve order in their map
        seqsToRemove = []
        for seq, block in blocks.items():
            if not self.isOrdered(block):
                f.write("seq %s is out of order\n" % seq)
            if self.isDuplicated(block):
                f.write("seq %s has self alignment\n" % seq)
            if not self.isOrdered(block) or self.isDuplicated(block):
                f.write("removing seq %s\n" % seq)
                seqsToRemove.append(seq)
        for seq in seqsToRemove:
            del blocks[seq]

        # take only the blocks from the target sequence with the most mapped
        # bases
        tSeqMapped = []
        for seq, value in blocks.items():
            tBlocks = map(itemgetter(1), value)
            mappedBlocks = reduce(lambda r, v: r + (v[1] - v[0]), tBlocks, 0)
            tSeqMapped.append((seq, mappedBlocks))
        if len(tSeqMapped) == 0:
            # can happen if the sequence doesn't map to the target at all
            return (None, None)
        tSeqName = max(tSeqMapped, key=itemgetter(1))[0]
        blocks = blocks[tSeqName]
        return (blocks, tStrands[tSeqName])

    def isOrdered(self, blocks):
        prev = None
        tBlocks = sorted(map(itemgetter(1), blocks), key=itemgetter(0))
        for block in tBlocks:
            if prev is not None:
                # Blocks are always in increasing order (if they
                # preserve order) since the query strand is always
                # positive and PSL blocks follow query order
                if prev[1] > block[0]:
                    return False
            prev = block
        return True

    def isDuplicated(self, blocks):
        """Checks query blocks for duplications"""
        blocks.sort(key=itemgetter(0))
        merged = []
        prev = None
        for (qBlock, _) in blocks:
            if prev is not None:
                if prev[1] > qBlock[0]:
                    return True
            prev = qBlock
        return False

    def isContiguousInTarget(self, bedLine):
        (blocks, tStrand) = self.liftover(bedLine)
        if blocks is None or tStrand is None:
            return False
        bedFields = bedLine.split()
        bedStart = int(bedFields[1])
        bedEnd = int(bedFields[2])

        bedIntrons = []
        if len(bedFields) == 12:
            blockStarts = map(int, bedFields[11].split(","))
            blockSizes = map(int, bedFields[10].split(","))
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

        qGaps = []
        tGaps = []
        prevqEnd = bedStart
        prevtEnd = None
        f = open('log', 'a')
        for (qBlock, tBlock) in blocks:
            qGap = qBlock[0] - prevqEnd
            tGap = 0
            if prevtEnd is not None:
                tGap = tBlock[0] - prevtEnd if tStrand == '+' else prevtEnd - tBlock[1]
            # Ignore any overlap of bed12 gaps (introns) and q/tGaps
            isIntron = False
            for intron in bedIntrons:
                f.write("%d %d\n" % (intron[0], intron[1]))
                # Bit hacky since this will still apply during exon skipping etc.
                if qBlock[0] >= intron[1] and prevqEnd <= intron[0]:
                    if tGap < qGap - self.maxIntronDiff or tGap > qGap + self.maxIntronDiff:
                        f.write("bad intron: %d %d\n" % (qGap, tGap))
                        return False
                    else:
                        qGaps.append(qGap - (intron[1] - intron[0]))
                        isIntron = True
            if not isIntron:
                qGaps.append(qGap)
                tGaps.append(tGap)
            prevqEnd = qBlock[1]
            prevtEnd = tBlock[1] if tStrand == '+' else tBlock[0]

        f.write("%s\n" % qGaps)
        f.write("%s\n" % tGaps)

        # Add up blocks and see if they are the required fraction of
        # the query sequence
        totalInBed = bedEnd - bedStart
        if len(bedFields) == 12:
            for intron in bedIntrons:
                totalInBed -= intron[1] - intron[0]
        total = 0
        for (qBlock, tBlock) in blocks:
            total += qBlock[1] - qBlock[0]
        if float(total)/totalInBed < self.requiredMapFraction:
            return False

        if self.noDeletions and len([i for i in qGaps if i > self.maxGap]) > 0:
            return False

        if len([i for i in tGaps if i > self.maxGap]) > 0:
            return False

        return True

    def getContiguousLines(self, bedPath, startLineNum=0, endLineNum=-1):
        for lineNum, line in enumerate(open(bedPath)):
            if lineNum < startLineNum:
                continue
            elif lineNum >= endLineNum:
                break

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
                        default=100, type=int)
    parser.add_argument("--noDeletions", help="care about deletion gaps",
                        default=False, action='store_true')
    parser.add_argument("--sliceNum", help="number of slices to create",
                        type=int, default=1)
    parser.add_argument("--maxIntronDiff", help="Maximum amount that intron "
                        "gaps are allowed to change by", type=int,
                        default=10000)
    parser.add_argument("--requiredMapFraction", help="Fraction of bases in "
                        "the query that need to map to the target to be "
                        "accepted", type=float, default=0.0)
    parser.add_argument("--printNumBases", help="instead of printing the "
                        "passing BED lines, print the number of bases that "
                        "passed if the line as a whole passed",
                        action='store_true', default=False)
    # TODO: option to allow dupes in the target
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
