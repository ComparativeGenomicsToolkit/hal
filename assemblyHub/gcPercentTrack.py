#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""Creating GC percent track for the hubs
"""
import os
from sonLib.bioio import system  
from jobTree.scriptTree.target import Target

class GetGCpercent( Target ):
    def __init__(self, genomedir, genome):
        Target.__init__(self)
        self.genomedir = genomedir
        self.genome = genome

    def run(self):
        twobitfile = os.path.join(self.genomedir, "%s.2bit" %self.genome)
        tempfile = os.path.join(self.genomedir, "%s.gc.wigVarStep.gz" %self.genome)
        cmd = "hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 %s %s | gzip -c > %s" %(self.genome, twobitfile, tempfile)
        system(cmd)
        chrsizefile = os.path.join(self.genomedir, "chrom.sizes")
        gcfile = os.path.join(self.genomedir, "%s.gc.bw" %self.genome)
        cmd = "wigToBigWig %s %s %s" %(tempfile, chrsizefile, gcfile)
        system(cmd)
        system("rm -f %s" %tempfile)

def writeTrackDb_gcPercent(f, genome):
    f.write("track gcPercent\n")
    f.write("longLabel GC Percent in 5-base Window\n")
    f.write("shortLabel GC Percent\n")
    f.write("type bigWig 0 100\n")
    f.write("group map\n")
    f.write("visibility dense\n")
    f.write("windowingFunction Mean\n")
    f.write("bigDataUrl %s.gc.bw\n" %genome)
    
    f.write("priority 2\n")
    f.write("autoScale Off\n")
    f.write("maxHeightPixels 128:36:16\n")
    f.write("graphTypeDefault Bar\n")
    f.write("gridDefault OFF\n")
    f.write("color 0,0,0\n")
    f.write("altColor 128,128,128\n")
    f.write("viewLimits 30:70\n")
    f.write("html ../documentation/gcPercent\n")
    f.write("\n")

def addGcOptions(parser):
    from optparse import OptionGroup
    #group = OptionGroup(parser, "GC PERCENT", "GC Percent in 5-base Window.")
    group = OptionGroup(parser, "GC PERCENT")
    group.add_option('--gcContent', dest='gcContent', action='store_true', default=False, help='If specified, make GC-content tracks. Default=%default')
    parser.add_option_group(group)

