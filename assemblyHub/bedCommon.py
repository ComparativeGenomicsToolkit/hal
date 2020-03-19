#!/usr/bin/env python3

'''
Mon Oct  7 15:27:23 PDT 2013
If a gene contains large intron, break it into separate bed entries
'''
import sys
from argparse import ArgumentParser
from sonLib.bioio import system  

#========= FILTER LONG INTRONS ==========
class BedFormatError(Exception):
    pass

class Bed():
    '''Bed record
    '''
    def __init__(self, line, tab, ucscNames):
        if tab:
            items = line.strip().split('\t')
        else:
            items = line.strip().split()
        if len(items) < 3: 
            raise BedFormatError("Bed format for this program requires a minimum of 3 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        self.tab = tab
        if ucscNames:
            self.chrom = items[0].split('.')[-1]
        else:
            self.chrom = items[0]
        self.chromStart = int(items[1]) #base 0
        self.chromEnd = int(items[2]) #exclusive
        if len(items) >= 4:
            self.name = items[3]
        self.bed12 = False

        if len(items) >= 12:
            self.bed12 = True
            self.score = items[4]
            self.strand = items[5]
            self.thickStart = int(items[6])
            self.thickEnd = int(items[7])
            self.itemRgb = items[8]
            self.blockCount = int(items[9])
            self.blockSizes = [ int(i) for i in items[10].rstrip(',').split(',') ]
            self.blockStarts = [ int(i) for i in items[11].rstrip(',').split(',') ]
            self.extra = items[12:]

    def __cmp__(self, other):
        if self.chrom != other.chrom:
            return cmp(self.chrom, other.chrom)
        elif self.chromStart != other.chromStart:
            return cmp(self.chromStart, other.chromStart)
        else:
            return cmp(self.chromEnd, other.chromEnd)

    def getStr12(self):
        s = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s" %\
               (self.chrom, self.chromStart, self.chromEnd, self.name, \
                self.score, self.strand, self.thickStart, self.thickEnd, \
                self.itemRgb, self.blockCount, ",".join([str(s) for s in self.blockSizes]), \
                ",".join([str(s) for s in self.blockStarts]))
        if len(self.extra) > 0:
            s += "\t" + "\t".join(self.extra)
        
        if not self.tab:
            s = s.replace("\t", " ")
        return s

def readBedFile(file, tab, ucscNames=True):
    beds = []
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        bed = Bed(line, tab, ucscNames)
        beds.append(bed)
    f.close()
    return beds

def splitBed(bed, i):
    import copy
    leftBed = copy.deepcopy(bed)
    rightBed = copy.deepcopy(bed)
    
    leftBed.chromEnd = bed.chromStart + bed.blockStarts[i] + bed.blockSizes[i]
    leftBed.thickStart = leftBed.chromStart
    leftBed.thickEnd = leftBed.chromEnd
    leftBed.blockCount = i + 1
    leftBed.blockStarts = bed.blockStarts[0:i+1]
    leftBed.blockSizes = bed.blockSizes[0:i+1]

    rightBed.chromStart = bed.chromStart + bed.blockStarts[i+1]
    rightBed.thickStart = rightBed.chromStart
    rightBed.thickEnd = rightBed.chromEnd
    rightBed.blockCount = bed.blockCount - (i + 1)
    rightBed.blockStarts = [s - bed.blockStarts[i+1] for s in bed.blockStarts[i+1:]]
    rightBed.blockSizes = bed.blockSizes[i+1:]
    return leftBed, rightBed

def filterLongIntrons_bed(bed, maxIntron):
    beds = []
    for i in range(0, bed.blockCount -1):
        intronStart = bed.blockStarts[i] + bed.blockSizes[i]
        intronEnd = bed.blockStarts[i+1]
        if intronEnd - intronStart > maxIntron:
            leftBed, rightBed = splitBed(bed, i)
            rightbeds = filterLongIntrons_bed(rightBed, maxIntron)
            beds.append(leftBed)
            beds.extend(rightbeds)
            break
    if len(beds) == 0:
        beds = [bed]
    return beds

def writeBeds12(f, beds):
    for b in beds:
        f.write("%s\n" %b.getStr12())

def filterLongIntrons(infile, outfile, maxIntron, tab, ucscNames=True):
    #If a "gene" contains long intron(s), break the bed entry into separate entries
    beds = readBedFile(infile, tab, ucscNames)
    if len(beds) == 0 or not beds[0].bed12:
        system("cp %s %s" %(infile, outfile))
    else:
        f = open(outfile, 'w')
        for bed in beds:
            newbeds = filterLongIntrons_bed(bed, maxIntron)
            writeBeds12(f, newbeds)
        f.close()

def tabifyBed(bedPath):
    """Overwrites the given space-separated bed file with a tab-separated
    version."""
    lines = []
    for line in open(bedPath):
        lines.append("\t".join(line.split(" ")))
    bedFile = open(bedPath, 'w')
    for line in lines:
        bedFile.write(line)

def untabifyBed(bedPath):
    """Overwrites the given tab-separated bed file with a space-separated
    version."""
    lines = []
    for line in open(bedPath):
        lines.append(" ".join(line.split("\t")))
    bedFile = open(bedPath, 'w')
    for line in lines:
        bedFile.write(line)

#def main():
#    usage = "%prog <in.bed> <out.bed> maxIntron"
#    parser = OptionParser(usage=usage)
#    options, args = parser.parse_args()
#    if len(args) < 3:
#        parser.error('Required 3 arguments. Only see %d.\n' %len(args))
#    infile = args[0]
#    outfile = args[1]
#    maxIntron = int(args[2])
#    filterLongIntrons(infile, outfile, maxIntron)
#
#if __name__ == '__main__':
#    main()
