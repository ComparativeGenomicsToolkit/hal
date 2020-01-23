#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen (nknguyen@soe.ucsc.edu)
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txtimport unittest

'''
Check liftover status
Input: bed_file (1 entry per element, e.g 1 line per gene)
       hal_file
       query_name
       target_name
       out_file
Output: print to out_file tab separated fields checking how each line
        of bed_file map to the target, including fields:
        <length>: number of bases of the original region in the bed_file
                  e.g for gene, it's the # non-intronic bases
        <mapped>: proportion of length that mapped to the target genome
        <rearrangment>: yes/no
        <in_frame>: yes/no
        <insertions>: list out the insertions
        <deletions>: list out the deletions
'''


import os
import sys
from optparse import OptionParser
from sets import Set

from sonLib.bioio import system
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack


class Status():
    def __init__(self, name):
        self.name = name
        self.length = -1
        self.map = 0
        self.ins = []
        self.dels = []
        self.oo = False 
        self.inframe = True

class Psl():
    '''Psl record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) != 21: 
            raise ValueError("Psl format requires 21 fields, line \n%s\n only has %d fields.\n" %(line, len(items)))
        self.desc = line
        self.matches = int(items[0])
        self.misMatches = int(items[1])
        self.repMatches = int(items[2])
        self.nCount = int(items[3])
        self.qNumInsert = int(items[4])
        self.qBaseInsert = int(items[5])  # number of bases inserted in query
        self.tNumInsert = int(items[6])  # number of inserts in target
        self.tBaseInsert = int(items[7]) #number of bases inserted in target
        self.strand = items[8]  # query strand
        self.qName = items[9] 
        self.qSize = int(items[10])
        self.qStart = int(items[11]) #base 0
        self.qEnd = int(items[12])
        self.tName = items[13]
        self.tSize = int(items[14])
        self.tStart = int(items[15])
        self.tEnd = int(items[16])
        self.blockCount = int(items[17])
        self.blockSizes = [int(s) for s in items[18].rstrip(',').split(',')]
        self.qStarts = [int(s) for s in items[19].rstrip(',').split(',')]
        self.tStarts = [int(s) for s in items[20].rstrip(',').split(',')]
        if len(self.blockSizes) != self.blockCount or len(self.qStarts) != self.blockCount or len(self.tStarts) != self.blockCount:
            raise ValueError("Psl format requires that the number of items in blockSizes, qStarts, tStarts is equal to blockCount. Line: %s\n" %line)

    def __cmp__(self, other):  # compare by query coordinate
        if self.qName != other.qName:
            return cmp(self.qName, other.qName)
        elif self.qStart != other.qStart:
            return cmp(self.qStart, other.qStart)
        else:
            return cmp(self.qEnd, other.qEnd)

class Bed():
    '''Bed record
    '''
    def __init__(self, line):
        items = line.strip().split('\t')
        if len(items) < 3: 
            raise BedFormatError(("Bed format for this program requires a " +
                   "minimum of 3 fields, line \n%s\n only has %d fields.\n" %
                   (line, len(items))))
        #self.chr = items[0]
        self.chr = items[0].split('.')[-1]
        try:
            self.start = int(items[1])  # base 0
            self.end = int(items[2])  # exclusive
        except ValueError:
            print("BED %s has wrong format\n" % line)
        self.name = ''
        if len(items) > 3:
            self.name = items[3]

        if len(items) >= 12:
            self.score = items[4]
            self.strand = items[5]
            assert self.strand == '-' or self.strand == '+'
            self.thickStart = int(items[6])  # base 0
            self.thickEnd = int(items[7])
            self.itemRgb = items[8]
            self.blockCount = int(items[9])
            self.blockSizes = [int(i) for i in items[10].rstrip(',').split(',')]
            self.blockStarts = [int(i) for i in items[11].rstrip(',').split(',')]
            assert len(self.blockSizes) == self.blockCount
            assert len(self.blockStarts) == self.blockCount
            
            #if blockStarts[0] != 0, convert start & end so that blockStarts[0] = 0
            if (len(self.blockStarts) > 0 and
                (self.blockStarts[0] != 0 or
                 self.end != self.start + self.blockStarts[-1] + self.blockSizes[-1])):
                offset = self.blockStarts[0]
                self.start += offset
                self.blockStarts = [s - offset for s in self.blockStarts]
                self.end = self.start + self.blockStarts[-1] + self.blockSizes[-1]
        else:
            self.score = '.'
            self.strand = '.'
            self.thickStart = self.start
            self.thickEnd = self.end
            self.itemRgb = '.'
            if len(items) >= 10:
                self.blockCount = int(items[9])
            else:
                self.blockCount = 1
            self.blockSizes = [self.end - self.start]
            self.blockStarts = [0]
        
    def __cmp__(self, other):
        if self.chr != other.chr:
            return cmp(self.chr, other.chr)
        elif self.start != other.start:
            return cmp(self.start, other.start)
        else:
            return cmp(self.end, other.end)

    def getStr(self):
        blockSizes = ','.join([str(s) for s in self.blockSizes])
        blockStarts = ','.join([str(s) for s in self.blockStarts])
        return "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s" \
               %(self.chr, self.start, self.end, self.name, self.score,\
                 self.strand, self.thickStart, self.thickEnd, self.itemRgb,\
                 self.blockCount, blockSizes, blockStarts)

def get_bed(file):
    f = open(file, 'r')
    lines = f.readlines()
    assert len(lines) == 1
    bed = Bed(lines[0])
    f.close()
    return bed

def psl_pos_target(psl):
    # make sure the target is on the positive strand
    if len(psl.strand) != 2 or psl.strand[1] != '-':
        return psl
    rvstrand = {'-': '+', '+': '-'}
    psl.strand = rvstrand[psl.strand[0]] + rvstrand[psl.strand[1]]
    sizes = []
    qstarts = []
    tstarts = []
    for i in range(psl.blockCount - 1, -1, -1):
        size = psl.blockSizes[i]
        qs = psl.qSize - (psl.qStarts[i] + size)
        ts = psl.tSize - (psl.tStarts[i] + size)
        sizes.append(size)
        qstarts.append(qs)
        tstarts.append(ts)
    psl.blockSizes = sizes
    psl.qStarts = qstarts
    psl.tStarts = tstarts
    return psl

def get_psls(file):
    psls = []
    f = open(file, 'r')
    for line in f:
        psl = Psl(line)
        psl = psl_pos_target(psl)
        psls.append(psl)
    f.close()
    return psls

def psl_getPosCoords(psl):
    # This function reverse the psl if query strand is -, to make sure that query strand is always +
    # make sure that the target strand of INPUT psl is always on the + strand
    assert len(psl.strand) < 2 or psl.strand[1] != '-'
    strand = psl.strand
    if psl.strand[0] == '-':
        qstarts = []
        tstarts = []
        sizes = []
        for i in range(psl.blockCount - 1, -1, -1):
            qstart = psl.qSize - (psl.qStarts[i] + psl.blockSizes[i])
            tstart = psl.tSize - (psl.tStarts[i] + psl.blockSizes[i])
            qstarts.append(qstart)
            tstarts.append(tstart)
            sizes.append(psl.blockSizes[i])
        qstrand = '+'
        if len(psl.strand) == 2 and psl.strand[1] == '-':
            tstrand = '+'
        else:
            tstrand = '-'
        strand = qstrand + tstrand
    else:
        qstarts = psl.qStarts
        tstarts = psl.tStarts
        sizes = psl.blockSizes
    return qstarts, tstarts, sizes, strand
    
def psl_check_query_overlap(psl1, psl2):
    # return True if query ranges of psl1 is overlap with query ranges of psl2
    overlap = 0
    if (psl1.qName != psl2.qName or psl1.qEnd <= psl2.qStart or
                                    psl2.qEnd <= psl1.qStart):  # not overlap
        return overlap 

    # Convert query coordinates of both psls to query + strand
    starts1, tstarts1, sizes1, strand1 = psl_getPosCoords(psl1)
    starts2, tstarts2, sizes2, strand2 = psl_getPosCoords(psl2)
    # Check each block:
    for i1, start1 in enumerate(starts1):
        end1 = start1 + sizes1[i1]
        for i2, start2 in enumerate(starts2):
            end2 = start2 + sizes2[i2]
            if start2 < end1 and start1 < end2:  # overlap
                ostart = max(start1, start2)
                oend = min(end1, end2)
                overlap += (oend - ostart)
    return overlap

def get_next_non_overlap_psl(sets, sorted_psls):
    new_sets = []  # list of lists, each elment = ([], lastindex)
    # each set is represented as list of indices of the psls in that list
    for curr_set in sets:
        psl_indices = curr_set[0]
        i = curr_set[1]  # current index
        added = 0
        for j in range(i + 1, len(sorted_psls)):
            psl = sorted_psls[j]
            overlap = False
            for index in psl_indices:
                psl0 = sorted_psls[index]
                if psl_check_query_overlap(psl, psl0) > 0:
                    overlap = True
                    break
            if not overlap:  # current psl does not overlap with any element
                added += 1
                new_set = (psl_indices + [j], j)
                curr_new_sets = get_next_non_overlap_psl([new_set], sorted_psls)
                new_sets = new_sets + curr_new_sets
        if added == 0:  # no additional non-overlap psl found
            new_set = (psl_indices, len(sorted_psls))
            new_sets.append(new_set)
    return new_sets 

def get_non_overlap_psls_sets(psls):
    # some bases in the query bed map to multiple positions on the target
    # this function return all possible sets of psls where the query base
    # only appear at most once
    sets = []
    for i in range(len(psls)):  # each psl as a starting element
        start_set = ([i], i)
        curr_sets = get_next_non_overlap_psl([start_set], psls)
        for s in curr_sets:  # make sure added set is not a subset of existing ones
            curr_set = Set(s[0])
            issubset = False
            for s0 in sets:
                set0 = Set(s0)
                if curr_set.issubset(set0):
                    issubset = True
                    break
            if not issubset:
                sets.append(s[0])
    return sets

def psls_get_qcov(psls):
    return sum([sum(psl.blockSizes) for psl in psls])

def get_most_qcov(pslsets, sorted_psls):
    # return the psl set with highest query coverage
    if not pslsets:
        return None
    most_qcov_set = None
    qcov = 0
    for pslset in pslsets:
        currcov = psls_get_qcov([sorted_psls[i] for i in pslset])
        if currcov > qcov:
            qcov = currcov
            most_qcov_set = pslset
    return most_qcov_set, qcov

class Reg:
    def __init__(self, name, start, end, strand, size, qstart, qend):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.size = size
        self.qstart = qstart
        self.qend = qend
        
def psl_get_tpos(qstarts, tstarts, sizes, qpos):
    for i, qs in enumerate(qstarts):
        size = sizes[i]
        qe = qs + size
        if qs <= qpos and qpos <= qe:
            offset = qpos - qs
            ts = tstarts[i]
            return ts + offset
    return -1

def block_status(indices, sorted_psls, start, end, edge):
    ins = []
    dels = []
    oo = False
    tregs = []
    tstart = -1
    tend = -1
    tname = ''

    pos = start
    blocksize = end - start
    assert blocksize > 0
    
    for index in indices:
        psl = sorted_psls[index]
        qstarts, tstarts, sizes, strand = psl_getPosCoords(psl)

        for i, qstart in enumerate(qstarts):
            qend = qstart + sizes[i]
            if qend < pos:
                continue
            elif end < qstart:
                break
            else:
                oqstart = max(pos, qstart)
                oqend = min(end, qend)
                otstart = psl_get_tpos(qstarts, tstarts, sizes, oqstart) 
                otend = psl_get_tpos(qstarts, tstarts, sizes, oqend)
                if strand[1] == '-':
                    temp = otstart
                    otstart = psl.tSize - otend
                    otend = psl.tSize - temp

                assert otend >= otstart
                #if otend < otstart:
                #    temp = otstart
                #    otstart = otend
                #    otend = otstart
                    
                treg = Reg(psl.tName, otstart, otend, strand[1], psl.tSize, oqstart, oqend)
                tregs.append(treg)
                
                if float(oqstart - start)/blocksize > edge:
                    deletion = oqstart - pos
                    if deletion > 0:
                        dels.append(deletion)
                pos = oqend
    if float(end - pos)/blocksize > edge:
        if pos < end:
            dels.append(end - pos)

    # checking for insertions:
    if len(tregs) > 1:
        for i in range(1, len(tregs)):
            if float(treg.qstart - start)/blocksize <= edge or float(end - treg.qend)/blocksize <= edge:
                continue
            prev_treg = tregs[i - 1]
            treg = tregs[i]
            if treg.name == prev_treg.name:
                if treg.strand == prev_treg.strand:
                    if treg.strand == '+':
                        if prev_treg.end < treg.start:
                            insertion = treg.start - prev_treg.end
                            ins.append(insertion)
                        elif prev_treg.end > treg.start:
                            oo = True
                    else:
                        if treg.end < prev_treg.start:
                            insertion = prev_treg.start - treg.end
                            ins.append(insertion)
                        elif treg.end > prev_treg.start:
                            oo = True
                else:
                    oo = True
            else:  # map to different chromosome
                oo = True
    
    strands = [treg.strand for treg in tregs]
    if len(tregs) > 0:
        tstart = min([treg.start for treg in tregs])
        tend = max([treg.end for treg in tregs])
        tname = tregs[0].name
    return ins, dels, oo, strands, tstart, tend, tname

def flipbed(bed):
    bed.blockSizes.reverse()
    bed.blockStarts.reverse()
    return bed

def get_liftover_status(bedfile, liftfile, edge):
    # return a status object
    bed = get_bed(bedfile)
    
    psls = get_psls(liftfile)
    status = Status(bed.name)
    # get length:
    l = sum(bed.blockSizes)
    status.length = l

    # get all possible mapping scenario (all sets of psls where each query base is unique)
    sorted_psls = sorted(psls)
    pslsets = get_non_overlap_psls_sets(sorted_psls)
    if not pslsets:
        return status

    most_qcov_set, map = get_most_qcov(pslsets, sorted_psls)

    # map, insertions, deletions, oo, inframe
    status.map = map
    ins = []
    dels = []
    currstrand = ''
    currtstart = -1
    currtend = -1
    currtname = ''
    for i, start in enumerate(bed.blockStarts):  # each block
        qstart = bed.start + start
        qend = qstart + bed.blockSizes[i]
        block_ins, block_dels, block_oo, strands, tstart, tend, tname = block_status(most_qcov_set, sorted_psls, qstart, qend, edge)
        ins.extend(block_ins)
        dels.extend(block_dels)
        if block_oo:
            status.oo = True
        elif strands: # check with previous block
            tstrand = strands[0]
            if currstrand:  # not first block
                if currstrand != tstrand and not status.oo:  # change orientation
                    status.oo = True
                elif currtname and tname and tname != currtname and not status.oo:
                    status.oo = True
                else:  #check for change in blocks' order
                    if ((tstrand == '+' and currtend > tstart) or
                        (tstrand == '-' and currtstart < tend)):
                        status.oo = True
            currstrand = tstrand
        if tstart > -1 and tend > -1:
            currtstart = tstart
            currtend = tend
            currtname = tname

    status.ins = ins
    status.dels = dels

    if status.oo or abs(sum(ins) - sum(dels)) % 3 > 0:
            status.inframe = False
    return status

def print_status(status, outfile):
    ins = ",".join([str(i) for i in status.ins])
    dels = ",".join([str(i) for i in status.dels])
    f = open(outfile, 'w')
    if status.map > 0:
        f.write("%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (status.name, status.length,
                 status.map, ins, dels, str(status.oo), str(status.inframe)))
    else:
        f.write("%s\t%d\t%d\t%s\t%s\tNA\tNA\n" % (status.name, status.length,
                                                  status.map, ins, dels))
    f.close()

def splitfile(file, outdir):
    f = open(file, 'r')
    i = 0
    for line in f:
        i += 1
        outfile = os.path.join(outdir, "%d.bed" % i)
        ofh = open(outfile, 'w')
        ofh.write(line)
        ofh.close()
    f.close()

class Setup(Target):
    def __init__(self, options):
        Target.__init__(self)
        self.opts = options
    
    def run(self):
        #Split bed file into separate entries
        global_dir = self.getGlobalTempDir()
        beddir = os.path.join(global_dir, "single_beds")
        system("mkdir -p %s" % beddir)
        splitfile(self.opts.bedfile, beddir)

        #For each bed, lift-over to the target genome, and report status
        liftoverdir = os.path.join(global_dir, "liftover")
        system("mkdir -p %s" % liftoverdir)
        outdir = os.path.join(global_dir, "out")
        system("mkdir -p %s" % outdir)
        for bedfile in os.listdir(beddir):
            bedpath = os.path.join(beddir, bedfile)
            liftoverfile = os.path.join(liftoverdir, bedfile)
            outfile = os.path.join(outdir, os.path.splitext(bedfile)[0])
            self.addChildTarget(LiftoverAndStatus(bedpath, liftoverfile,
                                                  outfile, self.opts))
        self.setFollowOnTarget(PrintResults(outdir, self.opts.outfile))

class PrintResults(Target):
    def __init__(self, indir, outfile):
        Target.__init__(self)
        self.indir = indir
        self.outfile = outfile
    
    def run(self):
        f = open(self.outfile, 'w')
        f.write("#Name\tLength\tMap\tIns\tDels\tOO\tInframe\n")
        f.close()
        system("cat %s/* >> %s" % (self.indir, self.outfile))

class LiftoverAndStatus(Target):
    def __init__(self, bedfile, liftoverfile, statusfile, opts):
        Target.__init__(self)
        self.bedfile = bedfile
        self.liftfile = liftoverfile
        self.statusfile = statusfile
        self.opts = opts

    def run(self):
        cmd = "halLiftover --outPSL --tab %s %s %s %s %s" % (self.opts.halfile,
                self.opts.query, self.bedfile, self.opts.target, self.liftfile)
        system(cmd)
        #system("cp %s %s_liftoverpsl" % (self.liftfile, self.opts.outfile))
        status = get_liftover_status(self.bedfile, self.liftfile, self.opts.edge)
        print_status(status, self.statusfile)

def addOptions(parser):
    parser.add_option('--edge', dest='edge', default=0.0,
                      help='proportion of block at each edge that is allowed to have errors')

def checkOptions(parser, args, options):
    if len(args) < 5:
        parser.error("Need 5 input arguments, only %d was given." % len(args))
    options.bedfile = args[0]
    if not os.path.exists(options.bedfile):
        parser.error("Input bed file %s does not exist." % options.bedfile)
    options.halfile = args[1]
    if not os.path.exists(options.halfile):
        parser.error("Input hal file %s does not exist." % options.halfile)
    options.query = args[2]
    options.target = args[3]
    options.outfile = args[4]

def main():
    usage = "%prog <bed_file> <hal_file> <query_name> <target_name> <out_file>"
    parser = OptionParser(usage=usage)
    addOptions(parser)
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    i = Stack(Setup(options)).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" % i)

if __name__ == '__main__':
    from halLiftoverStatus import *
    main()

