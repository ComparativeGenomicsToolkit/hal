#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Creating the Conservation track for the hubs
"""
import os
from sonLib.bioio import system  
from toil.job import Job
from optparse import OptionGroup

class GetConservationFiles( Job ):
    '''Compute conservation bigWig for the input genome and lift it over to all other genomes
    '''
    def __init__(self, halfile, outdir, options):
        Job.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        #First create the model:
        #halPhyloPTrain.py ../../out.hal reference genesMin66Samples.bed ref.mod --numProc 16 --tree rootedTree.nw
        modfile = os.path.join(self.outdir, "%s.mod" % self.options.conservationGenomeName)
        cmd = "halPhyloPTrain.py %s %s %s %s --numProc %d" %(self.halfile, \
                self.options.conservationGenomeName, self.options.conservation, \
                modfile, self.options.conservationNumProc)
        if self.options.conservationTree:
            cmd += " --tree %s" %self.options.conservationTree
        system(cmd)

        self.addFollowOn( GetConservationFiles2(self.halfile, self.outdir, modfile, self.options.conservationNumProc) )

class GetConservationFiles2( Job ):
    '''Modify mod file to convert small branch lengths (change all the xxxe-1y to xxxe-08)
       Then do liftover
    '''
    def __init__(self, halfile, outdir, modfile, numproc):
        Job.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.modfile = modfile
        self.numproc = numproc

    def run(self, fileStore):
        newmodfile = "%s-modified" %self.modfile
        #modify small branch lengths (change all the xxxe-1y to xxxe-10)
        system("sed 's/e-1./e-08/g' %s > %s" %(self.modfile, newmodfile))
        #get conservation bigwig and liftover files:
        cmd = "halTreePhyloP.py %s %s %s --bigWig --numProc %d" %(self.halfile, newmodfile, self.outdir, self.numproc)
        system(cmd)

def writeTrackDb_conservation(f, genome, conservationDir):
    wigfile = os.path.join(conservationDir, "%s_phyloP.bw" %genome)
    #if os.path.exists(wigfile):
    if True:#HACK
        f.write("track conservation\n")
        f.write("longLabel Conservation\n")
        f.write("shortLabel Conservation\n")
        f.write("type bigWig -1 1\n")
        f.write("group map\n")
        f.write("visibility dense\n")
        f.write("windowingFunction Mean\n")
        f.write("bigDataUrl ../conservation/%s_phyloP.bw\n" %genome)
        
        f.write("priority 2\n")
        f.write("autoScale On\n")
        f.write("maxHeightPixels 128:36:16\n")
        f.write("graphTypeDefault Bar\n")
        f.write("gridDefault OFF\n")
        f.write("color 0,0,0\n")
        f.write("altColor 128,128,128\n")
        f.write("viewLimits -1:1\n")
        f.write("html ../documentation/conservation\n")
        f.write("\n")

def addConservationOptions(parser):
    group = parser.add_argument_group("CONSERVATION TRACKS", "Necessary information for computing conservation tracks") 
    group.add_argument('--conservation', dest='conservation', help='Bed file providing regions to calculate the conservation tracks.')
    group.add_argument('--conservationDir', dest='conservationDir', help='Optional. Directory contains conservation bigwigs if available. These bigwigs will be used. If this is not specified, the program will compute the conservation tracks.')
    group.add_argument('--conservationGenomeName', dest='conservationGenomeName', help='Name of the genome of the bed file provided in the "--conversation" option')
    group.add_argument('--conservationTree', dest='conservationTree', help='Optional. Newick tree for the conservation track')
    group.add_argument('--conservationNumProc', dest='conservationNumProc', type=int, default=1, help='Optional. Number of processors to run conservation')
    group = parser.add_argument_group(group)

def checkConservationOptions(parser, options):
    if options.conservationDir:
        if not os.path.exists(options.conservationDir) or not os.path.isdir(options.conservationDir):
            parser.error("Conservation directory %s does not exist.\n" %options.conservationDir)
        options.conservation = True
    elif options.conservation:
        if not os.path.exists(options.conservation):
            parser.error("Conservation bed file %s does not exist.\n" %options.conservation)
        if not options.conservationGenomeName:
            parser.error("Please specify --conservationGenomeName, this option is required to compute the conservation tracks.\n")






