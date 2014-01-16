#!/usr/bin/env python
# won't work in the general case right now
# (i.e. inputs need to be single copy & positive query strand)
import sys
import argparse
from sonLib import bioio
from operator import itemgetter
from collections import defaultdict

class ContiguousRegions:
    def __init__(self, alignment, srcGenome, destGenome, maxGap):
        self.alignment = alignment
        self.srcGenome = srcGenome
        self.destGenome = destGenome
        self.maxGap = maxGap

    def liftover(self, bedLine):
        tempSrc = bioio.getTempFile("ContiguousRegions.tempSrc.bed")
        tempDest = bioio.getTempFile("ContiguousRegions.tempDest.psl")
        open(tempSrc, 'w').write("%s\n" % bedLine)
        cmd = "halLiftover --outPSL %s %s %s %s %s" % (self.alignment,
                                                       self.srcGenome,
                                                       tempSrc,
                                                       self.destGenome,
                                                       tempDest)
        bioio.system(cmd)
        pslLines = open(tempDest).read().split("\n")
        pslLines = map(lambda x: x.split(), pslLines)
        # dict is to keep blocks separated by target sequence
        qBlocks = defaultdict(list)
        tBlocks = defaultdict(list)
        for pslLine in pslLines:
            if pslLine == []:
                continue
            qStrand = pslLine[8][0]
            assert(qStrand == '+')
            tStrand = pslLine[8][1]
            tName = pslLine[13]
            blockSizes = [int(i) for i in pslLine[18].split(",") if i != '']
            qStarts = [int(i) for i in pslLine[19].split(",") if i != '']
            tStarts = [int(i) for i in pslLine[20].split(",") if i != '']
            assert(len(blockSizes) == len(qStarts) and
                   len(qStarts) == len(tStarts))
            for blockLen, qStart, tStart in zip(blockSizes, qStarts, tStarts):
                qBlocks[tName].append((qStart, qStart + blockLen))
                tBlocks[tName].append((tStart, tStart + blockLen))

        # take only the blocks from the target sequence with the most mapped
        # bases
        tSeqMapped = []
        for seq, blocks in tBlocks.items():
            mappedBlocks = reduce(lambda r, v: r + (v[1] - v[0]), blocks, 0)
            tSeqMapped.append((seq, mappedBlocks))
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

    def getContiguousLines(self, bedPath):
        for line in open(bedPath):
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
    parser.add_argument("--maxGap", help="maximum gap size to accept", 
                        default=10, type=int)
    args = parser.parse_args()

    contiguousRegions = ContiguousRegions(args.alignment, args.srcGenome,
                                          args.destGenome, args.maxGap)
    
    for line in contiguousRegions.getContiguousLines(args.bedFile):
        print line
    return 0
    
if __name__ == '__main__':
    sys.exit(main())
