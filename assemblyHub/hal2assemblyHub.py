#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen (nknguyen@soe.ucsc.edu)
#
#Released under the MIT license, see LICENSE.txtimport unittest

#Wed Apr 10 15:30:53 PDT 2013
#
#Generates necessary files to make assembly hub
#Input: 1/ Output directory 
#       2/ hal file of the multiple alignment
#       3/ (Optional: directory containing annotated bed files. e.g : genes)
#Output:
#   outdir/
#       hub.txt
#       genomes.txt
#
#http://genomewiki.ucsc.edu/index.php/Browser_Track_Construction

import os, sys, re, time
from optparse import OptionParser

from sonLib.bioio import system  
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

from hal.assemblyHub.prepareLodFiles import *
from hal.assemblyHub.prepareHubFiles import *
from hal.assemblyHub.alignabilityTrack import *
from hal.assemblyHub.bedTrack import *
from hal.assemblyHub.conservationTrack import *
from hal.assemblyHub.gcPercentTrack import *
from hal.assemblyHub.rmskTrack import *
from hal.assemblyHub.snakeTrack import *

###################### MAIN PIPELINE #####################
class Setup( Target ):
    '''Setting up the pipeline
    '''
    def __init__(self, halfile, outdir, options):
        Target.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self):
        writeHubFile(self.outdir, self.options)
        annotations = []
        if self.options.beddirs:
            annotations.extend( [os.path.basename(item) for item in self.options.beddirs.split(',')] )
        if self.options.bbdirs:
            annotations.extend( [os.path.basename(item) for item in self.options.bbdirs.split(',')] )
        writeGroupFile(self.outdir, annotations)
        allgenomes = getGenomesFromHal(self.halfile)
        genomes = []
        if self.options.genomes:
            for g in self.options.genomes:
                if g in allgenomes:
                    genomes.append(g)
        else:
            genomes = allgenomes
        genome2seq2len = getGenomeSequences(self.halfile, genomes)
        #Get basic files (2bit, chrom.sizes) for each genome:
        for genome in genomes: 
            self.addChildTarget( GetBasicFiles(genome, genome2seq2len[genome], self.halfile, self.outdir, self.options) )
        
        self.setFollowOnTarget( MakeTracks(genomes, genome2seq2len, self.halfile, self.outdir, self.options) )

class GetBasicFiles( Target ):
    '''Get 2bit and chrom.sizes for each genome
    '''
    def __init__(self, genome, seq2len, halfile, outdir, options):
        Target.__init__(self)
        self.genome = genome
        self.seq2len = seq2len
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self):
        genomedir = os.path.join(self.outdir, self.genome)
        system("mkdir -p %s" % genomedir)
        makeTwoBitSeqFile(self.genome, self.halfile, genomedir) #genomedir/genome.2bit
        getChromSizes(self.halfile, self.seq2len, os.path.join(genomedir, "chrom.sizes")) #genomedir/chrom.sizes

class MakeTracks( Target ):
    def __init__(self, genomes, genome2seq2len, halfile, outdir, options):
        Target.__init__(self)
        self.genomes = genomes
        self.genome2seq2len = genome2seq2len
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self):
        for genome in self.genomes:
            genomedir = os.path.join(self.outdir, genome)
            if self.options.gcContent:
                self.addChildTarget( GetGCpercent(genomedir, genome) ) #genomedir/genome.gc.bw
            if self.options.alignability:
                self.addChildTarget( GetAlignability(genomedir, genome, self.halfile) )#genomedir/genome.alignability.bw
        
        #Compute conservation track:
        #if not self.options.conservationDir and self.options.conservation: #HACK
        if self.options.conservation:
            conservationDir = os.path.join(self.outdir, "conservation")
            if not self.options.conservationDir: 
                system("mkdir -p %s" %conservationDir)
                self.addChildTarget( GetConservationFiles(self.halfile, conservationDir, self.options) )
            else:
                if os.path.abspath(self.options.conservationDir) != os.path.abspath(conservationDir):
                    system("cp -r %s %s" %(self.options.conservationDir, conservationDir))

        #Lift-over annotations of each sample to all other samples
        self.options.bigbeddirs = []
        if self.options.beddirs:
            beddirs = self.options.beddirs.split(',')
            for beddir in beddirs:
                annotation = os.path.basename(beddir)
                bigbeddir = os.path.join(self.outdir, "liftoverbeds", annotation)
                system("mkdir -p %s" % bigbeddir)
                self.options.bigbeddirs.append(bigbeddir)
                self.addChildTarget( LiftoverBedFiles(beddir, self.halfile, self.genome2seq2len, bigbeddir, self.options.noBedLiftover, self.outdir) )
        #Add the 'final' annotation tracks, no liftover
        if self.options.bbdirs:
            bbdirs = self.options.bbdirs.split(',')
            liftoverdir = os.path.join(self.outdir, "liftoverbeds")
            if not os.path.exists(liftoverdir):
                system("mkdir -p %s" %liftoverdir)
            for bbdir in bbdirs:
                annotation = os.path.basename(bbdir)
                bigbeddir = os.path.join(self.outdir, "liftoverbeds", annotation)
                self.options.bigbeddirs.append(bigbeddir)
                #Copy bb files to bigbeddir
                if os.path.abspath(bigbeddir) != os.path.abspath(bbdir):
                    #system("mkdir -p %s" %bigbeddir)
                    system("cp -r %s %s" %(bbdir, bigbeddir))

        #Get LOD if needed, and Write trackDb files
        self.setFollowOnTarget( WriteGenomesFile(self.genomes, self.genome2seq2len, self.halfile, self.options, self.outdir) )

class WriteGenomesFile(Target):
    '''Write genome for all samples in hal file
    '''
    def __init__(self, genomes, genome2seq2len, halfile, options, outdir):
        Target.__init__(self)
        self.genomes = genomes
        self.genome2seq2len = genome2seq2len
        self.halfile = halfile
        self.options = options
        self.outdir = outdir

    def run(self):
        options = self.options
        localHalfile = os.path.join(self.outdir, os.path.basename(self.halfile))
        if os.path.abspath(localHalfile) != os.path.abspath(self.halfile):
            if os.path.exists(localHalfile):
                system("rm %s" %localHalfile)
            if options.cpHal:
                system("cp %s %s" %(os.path.abspath(self.halfile), localHalfile))
            else:
                system("ln -s %s %s" %(os.path.abspath(self.halfile), localHalfile))

        #Create lod files if useLod is specified
        lodtxtfile, loddir = getLod(options, localHalfile, self.outdir)
        
        filename = os.path.join(self.outdir, "genomes.txt")
        f = open(filename, 'w')
        for genome in self.genomes:
            genomedir = os.path.join(self.outdir, genome)
            f.write("genome %s\n" %genome)
            f.write("twoBitPath %s/%s.2bit\n" % (genome, genome))

            #create trackDb for the current genome:
            if lodtxtfile == '':
                self.addChildTarget( WriteTrackDbFile(self.genomes, "../%s" % os.path.basename(self.halfile), genomedir, options) )
            else:
                self.addChildTarget( WriteTrackDbFile(self.genomes, "../%s" % os.path.basename(lodtxtfile), genomedir, options) )
            f.write("trackDb %s/trackDb.txt\n" %genome)
            
            #other info
            f.write("groups groups.txt\n")

            writeDescriptionFile(genome, genomedir)
            f.write("htmlPath %s/description.html\n" %genome)
            f.write("organism %s\n" %genome)
            f.write("orderKey 4800\n")
            f.write("scientificName %s\n" %genome)
            
            seq2len = self.genome2seq2len[genome]
            (seq, l) = getLongestSeq(seq2len)
            f.write("defaultPos %s:1-%d\n" %(seq, min(l, 1000)))
            f.write("\n")
        f.close()

class WriteTrackDbFile( Target ):
    def __init__(self, genomes, halfile, outdir, options):
        Target.__init__(self)
        self.genomes = genomes
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self):
        currgenome = self.outdir.rstrip('/').split("/")[-1]
        filename = os.path.join(self.outdir, "trackDb.txt")
        f = open(filename, 'w')
        
        if self.options.gcContent:
            writeTrackDb_gcPercent(f, currgenome)
        if self.options.alignability:
            writeTrackDb_alignability(f, currgenome, len(self.genomes))
        if self.options.conservation:
            conservationDir = os.path.join(self.outdir, "..", "conservation")
            writeTrackDb_conservation(f, currgenome, conservationDir)

        for bigbeddir in self.options.bigbeddirs:
            writeTrackDb_bigbeds(f, bigbeddir, self.genomes, currgenome, self.options.properName)

        if self.options.rmskdir:
            writeTrackDb_rmsk(f, os.path.join(self.options.rmskdir, currgenome), self.outdir)

        writeTrackDb_snakes(f, self.halfile, self.genomes, currgenome, self.options.properName)
        f.close()

############################ UTILITIES FUNCTIONS ###################
def getLongestSeq(seq2len):
    seqs = sorted( [(seq, len) for seq, len in seq2len.iteritems()], key=lambda item:item[1], reverse=True )
    return seqs[0]

def getGenomeSequencesFromHal(halfile, genome):
    statsfile = "%s-seqStats.txt" %genome
    system("halStats --sequenceStats %s %s > %s" %(genome, halfile, statsfile))
    
    seq2len = {}
    f = open(statsfile, 'r')
    for line in f:
        if len(line) < 2 or re.search("SequenceName", line):
            continue
        items = line.strip().split(", ")
        seq = items[0].split('.')[-1]
        #seq = items[0]
        l = int(items[1])
        seq2len[seq] = l
    f.close()
    system("rm %s" %statsfile)

    return seq2len

def getGenomeSequences(halfile, genomes):
    genome2seq2len = {}
    for genome in genomes:
        seq2len = getGenomeSequencesFromHal(halfile, genome)
        if len(seq2len) == 0:
            sys.stderr.write("Warning: genome %s contains 0 sequence - no browser was made.\n" %genome)
        else:
            genome2seq2len[genome] = seq2len
    return genome2seq2len

def getChromSizes(halfile, seq2len, outfile):
    f = open(outfile, 'w')
    for s, l in seq2len.iteritems():
        if l > 0:
            f.write("%s\t%d\n" %(s, l))
    f.close()

def getGenomesFromHal(halfile):
    #Get a list of all genomes from the output of halStats
    statsfile = "halStats.txt"
    system("halStats --genomes %s > %s" %(halfile, statsfile))
    
    f = open(statsfile, 'r')
    genomes = f.readline().strip().split()
    f.close()

    #clean up
    system("rm %s" %statsfile)
    
    return genomes

def makeTwoBitSeqFile(genome, halfile, outdir):
    fafile = os.path.join(outdir, "%s.fa" %genome)
    system("hal2fasta --outFaPath %s %s %s" %(fafile, halfile, genome))
    
    #if sequence headers have "." (e.g genome.chr), reformat the header to only have "chr"
    fafile2 = "%s2" %fafile
    cmd = "awk '{ if($0 ~/>/){split($1, arr, \".\"); if(length(arr) > 1 ){print \">\" arr[2]}else{print $0} }else{ print $0} }' %s > %s" %(fafile, fafile2)
    system(cmd)
    system("rm %s" %fafile)

    #convert to 2bit files
    twobitfile = os.path.join(outdir, "%s.2bit" %genome)
    system("faToTwoBit %s %s" %(fafile2, twobitfile))
    system("rm %s" %fafile2)

def getFilesByExt(indir, fileExtension):
    #Return all files in indir that have "fileExtension"
    allfiles = os.listdir(indir)
    files = []
    for f in allfiles:
        if os.path.splitext(f) == ".%s" %fileExtension :
            files.append(f)
    return files

def addOptions(parser):
    parser.add_option('--cpHalFileToOut', dest='cpHal', action='store_true', default=False, help='If specified, copy the input halfile to the output directory (instead of just make a softlink). Default=%default')
    addHubOptions(parser)
    addLodOptions(parser)
    addBedOptions(parser)
    addRmskOptions(parser)
    addGcOptions(parser)
    addAlignabilityOptions(parser)
    addConservationOptions(parser)

def checkOptions(parser, args, options):
    if len(args) < 2:
        parser.error("Required two input arguments, %d was provided\n" %len(args))
    if not os.path.exists(args[0]):
        parser.error("Input hal file %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        system("mkdir -p %s" %args[1])
    elif not os.path.isdir(args[1]):
        parser.error("Output directory specified (%s) is not a directory\n" %args[1])
    checkHubOptions(parser, options)
    checkRmskOptions(parser, options)
    checkConservationOptions(parser, options)

def main():
    usage = '%prog <halFile> <outputDirectory> [options]'
    parser = OptionParser(usage = usage)
    addOptions(parser)
    Stack.addJobTreeOptions(parser)

    options, args = parser.parse_args()
    checkOptions(parser, args, options)
    
    halfile = args[0]
    outdir = args[1]

    i = Stack( Setup(halfile, outdir, options) ).startJobTree(options)
    if i:
        raise RuntimeError("The jobtree contains %d failed jobs.\n" %i)

if __name__ == '__main__':
    from hal.assemblyHub.hal2assemblyHub import *
    main()

