#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Creating Alignability (Alignment Depth) track for the hubs
"""
import os
from sonLib.bioio import system
from toil.job import Job

class GetAlignability( Job ):
    def __init__(self, genomedir, genome, halfile):
        Job.__init__(self)
        self.genomedir = genomedir
        self.genome = genome
        self.halfile = halfile

    def run(self, fileStore):
        outfile = os.path.join(self.genomedir, "%s.alignability.bw" %self.genome)
        tempwig = os.path.join(self.genomedir, "%s.alignability.wig" %self.genome)
        system("halAlignmentDepth %s %s > %s" %(self.halfile, self.genome, tempwig))
        chrsizefile = os.path.join(self.genomedir, "chrom.sizes")
        system("wigToBigWig %s %s %s" %(tempwig, chrsizefile, outfile))
        system("rm -f %s" %tempwig)

def writeTrackDb_alignability(f, genome, genomeCount):
    f.write("track alignability\n")
    f.write("longLabel Alignability\n")
    f.write("shortLabel Alignability\n")
    f.write("type bigWig 0 %d\n" %genomeCount)
    f.write("group map\n")
    f.write("visibility dense\n")
    f.write("windowingFunction Mean\n")
    f.write("bigDataUrl %s.alignability.bw\n" %genome)
    
    f.write("priority 2\n")
    f.write("autoScale On\n")
    f.write("maxHeightPixels 128:36:16\n")
    f.write("graphTypeDefault Bar\n")
    f.write("gridDefault OFF\n")
    f.write("color 0,0,0\n")
    f.write("altColor 128,128,128\n")
    f.write("viewLimits 0:%d\n" %genomeCount)
    f.write("html ../documentation/alignability\n")
    f.write("\n")

def addAlignabilityOptions(parser):
    from optparse import OptionGroup
    #group = parser.add_argument_group("ALIGNABILITY", "Alignability: the number of genomes aligned to each position.")
    group = parser.add_argument_group("ALIGNABILITY")
    group.add_argument('--alignability', dest='alignability', action='store_true', default=False, help='If specified, make Alignability (aka Alignment Depth) tracks. ')
    group = parser.add_argument_group(group)

