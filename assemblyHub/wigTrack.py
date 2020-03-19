#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Creating wiggle (annotation) tracks and lifted-over wiggle tracks for the hubs
"""
import os, re, time
from sonLib.bioio import system  
from toil.job import Job
from optparse import OptionGroup
from hal.assemblyHub.assemblyHubCommon import *

class LiftoverWigFiles( Job ):
    def __init__(self, indir, halfile, genome2seq2len, bigwigdir, noLiftover, outdir):
        Job.__init__(self)
        self.indir = indir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.bigwigdir = bigwigdir
        self.noLiftover = noLiftover
        self.outdir = outdir

    def run(self, fileStore):
        #wigdir has the hierachy: indir/genome/chr1.wig, chr2.wig...
        #for each genome in wigdir, liftover the wig records of that genome to the coordinate of all other genomes
         
        #liftover wig file of each genome with available wigs to all genomes
        genomes = list(self.genome2seq2len.keys())
        tempwigs = []
        
        for genome in os.listdir(self.indir):
            if genome not in genomes:
                continue
            genomeindir = os.path.join(self.indir, genome)
            assert os.path.isdir(genomeindir)

            #Create wig directory for current genome
            genomeoutdir = os.path.join(self.bigwigdir, genome)
            system("mkdir -p %s" %genomeoutdir)
        
            #get all the wig files (".wig" ext)
            wigfiles = getFilesByExt(genomeindir, "wig")

            #Concatenate all the input wig files and convert it into bigwig to outdir/genome/genome.bw
            tempwig = "%s-temp.wig" % os.path.join(genomeoutdir, genome)
            system( "cat %s/*wig > %s" %(genomeindir, tempwig) )
            if os.stat(tempwig).st_size > 0:#make sure the file is not empty
                outbigwig = os.path.join(genomeoutdir, "%s.bw" %genome)
                chrsizefile = os.path.join(self.outdir, genome, "chrom.sizes")
                system("wigToBigWig %s %s %s" %(tempwig, chrsizefile, outbigwig))
                
                #Liftover to all other genomes:
                if not self.noLiftover:
                    for othergenome in genomes:
                        if othergenome != genome:
                            self.addChild( LiftoverWig(genomeoutdir, tempwig, genome, othergenome, self.halfile, self.outdir) )
            tempwigs.append( tempwig )
        self.addFollowOn( CleanupFiles(tempwigs) )

class LiftoverWig( Job ):
    def __init__(self, genomeoutdir, wig, genome, othergenome, halfile, outdir):
        Job.__init__(self)
        self.genomeoutdir = genomeoutdir
        self.wig = wig
        self.genome = genome
        self.othergenome = othergenome
        self.halfile = halfile
        self.outdir = outdir

    def run(self, fileStore):
        liftovertempwig = "%s.wig" % os.path.join(self.genomeoutdir, self.othergenome)
        system("halWiggleLiftover %s %s %s %s %s" %(self.halfile, self.genome, self.wig, self.othergenome, liftovertempwig))
        outbigwig = os.path.join(self.genomeoutdir, "%s.bw" %self.othergenome)
        chrsizefile = os.path.join(self.outdir, self.othergenome, "chrom.sizes")
        if os.stat(liftovertempwig).st_size > 0:#make sure the file is not empty
            system("wigToBigWig %s %s %s" %(liftovertempwig, chrsizefile, outbigwig))
        #Cleanup:
        system("rm %s" % liftovertempwig)

#def writeTrackDb_bigwigs(f, bigwigdir, genomes, subgenomes, currgenome, properName):
def writeTrackDb_bigwigs(f, bigwigdir, genomes, currgenome, properName):
    annotation = os.path.basename(bigwigdir)
    genome2priority = {}
    for i, genome in enumerate(genomes):
        if genome == currgenome:
            genome2priority[genome] = 1
        else:
            genome2priority[genome] = i + 2
    
    for genome in os.listdir(bigwigdir):
        bwfile = os.path.join(bigwigdir, genome, "%s.bw" %currgenome)
        if not os.path.exists(bwfile):
            continue
        #start writing track
        genomeProperName = genome
        if genome in properName:
            genomeProperName = properName[genome]
        priority = 1
        if genome in genome2priority:
            priority = genome2priority[genome]

        f.write("\t\ttrack %s%s\n" % (annotation, genome))
        if genome == currgenome:
            f.write("\t\tlongLabel %s %s\n" % (genomeProperName, annotation))
        else:
            f.write("\t\tlongLabel %s Lifted-over %s\n" % (genomeProperName, annotation))
        f.write("\t\tpriority %d\n" %priority)
        f.write("\t\tshortLabel %s%s\n" % (genomeProperName, annotation))
        f.write("\t\tbigDataUrl ../liftoverwig/%s\n" % os.path.join( annotation, genome, "%s.bw" % currgenome ) )
        f.write("\t\ttype bigWig\n")
        f.write("\t\tgroup annotation%s\n" %annotation)
        f.write("\t\titemRgb On\n")
        #if genome == currgenome or genome in subgenomes:
        if genome == currgenome:
            f.write("\t\tvisibility dense\n")
            f.write("\t\tparent hubCentral%s\n"%annotation)
        else:
            f.write("\t\tvisibility hide\n")
            f.write("\t\tparent hubCentral%s off\n"%annotation)
        f.write("\t\twindowingFunction Mean\n")
        f.write("\t\tautoScale On\n")
        f.write("\t\tmaxHeightPixels 128:36:16\n")
        f.write("\t\tgraphTypeDefault Bar\n")
        f.write("\t\tgridDefault OFF\n")
        f.write("\t\tcolor 0,0,0\n")
        f.write("\t\taltColor 128,128,128\n")
        f.write("\t\tviewLimits 30:70\n")
        f.write("\t\tsubGroups view=%s orgs=%s\n" %(annotation, genome))
        f.write("\n")

def addWigOptions(parser):
    group = parser.add_argument_group("WIGGLE-FORMATTED ANNOTATIONS", "All annotations in wiggle or bigWig formats.")
    group.add_argument('--wigDirs', dest='wigdirs', help='comma separated list of directories containing wig files of the input genomes. Each directory represents a type of annotation. The annotations of each genome will then be liftovered to all other genomes in the MSA. Example: "genes,genomicIsland,tRNA". Format of each directory: wigDir/ then genome1/ then chr1.wig, chr2.wig... ' )
    group.add_argument('--finalBigwigDirs', dest='bwdirs', help='comma separated list of directories containing final big wig files to be displayed. No liftover will be done for these files. Each directory represents a type of annotation. Example: "readCoverage,". Format of each directory: bwDir/ then queryGenome/ then targetGenome1.bw, targetGenome2.bw ... (so annotation of queryGenome has been mapped to targetGenomes and will be display on the targetGenome browsers). ' )
    group.add_argument('--nowigLiftover', dest='noWigLiftover', action='store_true', default=False, help='If specified, will not lift over the wig annotations. ')
    group = parser.add_argument_group(group)

def checkWigOptions(parser, options):
    options.bigwigdirs = []
    if options.wigdirs:
        dirs = [d.rstrip('/') for d in options.wigdirs.split(',')]
        options.wigdirs = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Wig directory %s does not exist or is not a directory.\n" %d)
    if options.bwdirs:
        dirs = [d.rstrip('/') for d in options.bwdirs.split(',')]
        options.bwdirs = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Bigwig directory %s does not exist or is not a directory.\n" %d)




