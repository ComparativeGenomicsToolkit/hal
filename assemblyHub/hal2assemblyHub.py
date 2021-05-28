#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen (nknguyen@soe.ucsc.edu)
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
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
from argparse import ArgumentParser

from sonLib.bioio import system, getTempFile
from toil.job import Job
from toil.common import Toil

from hal.assemblyHub.prepareLodFiles import *
from hal.assemblyHub.prepareHubFiles import *
from hal.assemblyHub.alignabilityTrack import *
from hal.assemblyHub.bedTrack import *
from hal.assemblyHub.wigTrack import *
from hal.assemblyHub.conservationTrack import *
from hal.assemblyHub.gcPercentTrack import *
from hal.assemblyHub.groupExclusiveRegions import *
from hal.assemblyHub.rmskTrack import *
from hal.assemblyHub.snakeTrack import *
from hal.assemblyHub.assemblyHubCommon import MakeAnnotationTracks, preprocessAnnotationInputs, getProperName, sortByProperName
from hal.assemblyHub.treeCommon import *
from hal.assemblyHub.docs.makeDocs import *

###################### MAIN PIPELINE #####################
class Setup( Job ):
    '''Setting up the pipeline
    '''
    def __init__(self, halfile, outdir, options):
        Job.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        writeHubFile(self.outdir, self.options)
        annotations = []
        if self.options.beddirs:
            annotations.extend( [os.path.basename(item) for item in self.options.beddirs] )
        if self.options.bbdirs:
            annotations.extend( [os.path.basename(item) for item in self.options.bbdirs] )
        
        if self.options.beddirs2:
            annotations.extend( [os.path.basename(item) for item in self.options.beddirs2] )
        if self.options.bbdirs2:
            annotations.extend( [os.path.basename(item) for item in self.options.bbdirs2] )
        
        if self.options.wigdirs:
            annotations.extend( [os.path.basename(item) for item in self.options.wigdirs] )
        if self.options.bwdirs:
            annotations.extend( [os.path.basename(item) for item in self.options.bwdirs] )
        writeGroupFile(self.outdir, self.options.longLabel, annotations)
        # Check that no two annotations have the same name
        if not all([annotations.count(x) == 1 for x in annotations]):
            raise RuntimeError("Some annotation folders (beds, wigs) have the "
                               "same basename. Since the name for each "
                               "annotation track is taken from the basename "
                               "of the annotation folder, all basenames must "
                               "be unique.")

        #Get tree
        if not self.options.tree:
            checkHalTree(self.halfile, self.outdir, self.options)
        assert self.options.tree is not None # if this goes wrong yell at joel
        if isBinaryTree(self.options.tree): #get the png of the tree
            self.options.treeFig, self.options.leaves = drawTreeWtInternalNodesAligned(self.options.tree, self.outdir, self.options.properName)
        else:
            # Can't get tree png
            self.options.leaves = getLeaves(self.options.tree)

        #Get the ordering of the tracks
        #getOrderFromTree(self.options)
        allgenomes = getGenomesFromHal(self.halfile)
        genomes = []
        if self.options.genomes:
            for g in self.options.genomes:
                if g in allgenomes:
                    genomes.append(g)
        else:
            genomes = allgenomes
        genome2seq2len = getGenomeSequences(self.halfile, genomes,
                                            self.options.ucscNames)
        #Get basic files (2bit, chrom.sizes) for each genome:
        for genome in genomes: 
            self.addChild( GetBasicFiles(genome, genome2seq2len[genome], self.halfile, self.outdir, self.options) )
        
        self.addFollowOn( MakeTracks(genomes, genome2seq2len, self.halfile, self.outdir, self.options) )

class GetBasicFiles( Job ):
    '''Get 2bit and chrom.sizes for each genome
    '''
    def __init__(self, genome, seq2len, halfile, outdir, options):
        Job.__init__(self)
        self.genome = genome
        self.seq2len = seq2len
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        genomedir = os.path.join(self.outdir, self.genome)
        system("mkdir -p %s" % genomedir)
        if not self.options.twobitdir:
            self.makeTwoBitSeqFile(genomedir) #genomedir/genome.2bit
        else:
            linkTwoBitSeqFile(self.genome, self.options.twobitdir, genomedir) #genomedir/genome.2bit
        if self.options.ucscNames:
            # allows renaming of sequences
            getChromSizes(self.halfile, self.seq2len, os.path.join(genomedir, "chrom.sizes")) #genomedir/chrom.sizes
        else:
            getChromSizesFromHal(self.halfile, self.genome, os.path.join(genomedir, "chrom.sizes"))

    def makeTwoBitSeqFile(self, outdir):
        fafile = os.path.join(outdir, "%s.fa" %self.genome)
        system("hal2fasta --outFaPath %s %s %s" %(fafile, self.halfile, self.genome))
        fafile2 = fafile
        if self.options.ucscNames:
            #if sequence headers have "." (e.g genome.chr), reformat
            #the header to only have "chr"
            fafile2 = "%s2" %fafile
            cmd = "awk '{ if($0 ~/>/){split($1, arr, \".\"); if(length(arr) > 1 ){print \">\" arr[2]}else{print $0} }else{ print $0} }' %s > %s" %(fafile, fafile2)
            system(cmd)
            system("rm %s" %fafile)

        #convert to 2bit files
        twobitfile = os.path.join(outdir, "%s.2bit" %self.genome)
        system("faToTwoBit %s %s" %(fafile2, twobitfile))
        system("rm %s" %fafile2)

class MakeTracks( Job ):
    def __init__(self, genomes, genome2seq2len, halfile, outdir, options):
        Job.__init__(self)
        self.genomes = genomes
        self.genome2seq2len = genome2seq2len
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        #GC content & Alignability
        for genome in self.genomes:
            genomedir = os.path.join(self.outdir, genome)
            if self.options.gcContent:
                self.addChild( GetGCpercent(genomedir, genome) ) #genomedir/genome.gc.bw
            if self.options.alignability:
                self.addChild( GetAlignability(genomedir, genome, self.halfile) )#genomedir/genome.alignability.bw
        
        #Compute conservation track:
        if self.options.conservation:
        #if self.options.conservation or self.options.conservationDir:
            conservationDir = os.path.join(self.outdir, "conservation")
            if not self.options.conservationDir: 
                system("mkdir -p %s" %conservationDir)
                self.addChild( GetConservationFiles(self.halfile, conservationDir, self.options) )
            else:
                if os.path.abspath(self.options.conservationDir) != os.path.abspath(conservationDir):
                    system("ln -s %s %s" %(os.path.abspath(self.options.conservationDir), conservationDir))
                    #system("cp -r %s %s" %(self.options.conservationDir, conservationDir))

        #Make bed tracks:
        preprocessAnnotationInputs(self.options, self.outdir, "bed") 
        self.addChild( MakeAnnotationTracks(self.options, self.outdir, self.halfile, self.genome2seq2len, "bed") )
        
        #Make bed2 tracks:
        preprocessAnnotationInputs(self.options, self.outdir, "bed2") 
        self.addChild( MakeAnnotationTracks(self.options, self.outdir, self.halfile, self.genome2seq2len, "bed2") )
        
        #Make wig tracks:
        preprocessAnnotationInputs(self.options, self.outdir, "wig") 
        self.addChild( MakeAnnotationTracks(self.options, self.outdir, self.halfile, self.genome2seq2len, "wig") )

        #Make clade-exclusive tracks:
        if self.options.tree and self.options.cladeExclusive:
            self.addChild(GetCladeExclusiveRegions(self.halfile, self.options.tree, os.path.join(self.outdir, "liftoverbeds"), self.options.maxOut, self.options.minIn))
            self.options.bigbeddirs.append( os.path.join(self.outdir, "liftoverbeds", "CladeExclusive") )

        #Get LOD if needed, and Write trackDb files
        self.addFollowOn( WriteGenomesFile(self.genomes, self.genome2seq2len, self.halfile, self.options, self.outdir) )

class WriteGenomesFile(Job):
    '''Write genome for all samples in hal file
    '''
    def __init__(self, genomes, genome2seq2len, halfile, options, outdir):
        Job.__init__(self)
        self.genomes = genomes
        self.genome2seq2len = genome2seq2len
        self.halfile = halfile
        self.options = options
        self.outdir = outdir

    def run(self, fileStore):
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
        
        #Get the maximum window size to display SNPs
        if lodtxtfile:
            snpwidth = getLodLowestLevel(lodtxtfile) - 1
            if snpwidth > -1:
                options.snpwidth = snpwidth

        genomes = sortByProperName(self.genomes, self.options.properName)

        #Create documentation files:
        docdir = os.path.join(self.outdir, "documentation")
        system("mkdir -p %s" %docdir)
        writeDocFiles(docdir, self.options)

        #Create genomes.txt file
        filename = os.path.join(self.outdir, "genomes.txt")
        f = open(filename, 'w')
        #for genome in self.genomes:
        for genome in genomes:
            genomedir = os.path.join(self.outdir, genome)
            f.write("genome %s\n" %genome)
            f.write("twoBitPath %s/%s.2bit\n" % (genome, genome))

            #create trackDb for the current genome:
            if lodtxtfile == '':
                self.addChild( WriteTrackDbFile(self.genomes, "../%s" % os.path.basename(self.halfile), genomedir, options) )
            else:
                self.addChild( WriteTrackDbFile(self.genomes, "../%s" % os.path.basename(lodtxtfile), genomedir, options) )
            f.write("trackDb %s/trackDb.txt\n" %genome)
            
            #other info
            f.write("groups groups.txt\n")

            writeDescriptionFile(genome, genomedir)
            f.write("htmlPath %s/description.html\n" % genome)
            f.write("description %s\n" % getProperName(genome, self.options.properName))
            f.write("organism %s\n" % getProperName(genome, self.options.properName))
            f.write("orderKey 4800\n")
            f.write("scientificName %s\n" % genome)
            
            seq2len = self.genome2seq2len[genome]
            (seq, l) = getLongestSeq(seq2len)
            f.write("defaultPos %s:1-%d\n" %(seq, min(l, 1000)))
            f.write("\n")
        f.close()

class WriteTrackDbFile( Job ):
    def __init__(self, genomes, halfile, outdir, options):
        Job.__init__(self)
        self.genomes = genomes
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        currgenome = self.outdir.rstrip('/').split("/")[-1]
        filename = os.path.join(self.outdir, "trackDb.txt")
        f = open(filename, 'w')
        
        if self.options.gcContent:
            writeTrackDb_gcPercent(f, currgenome)
        if self.options.alignability:
            writeTrackDb_alignability(f, currgenome, len(self.genomes))
        if self.options.conservation:
        #if self.options.conservation or self.options.conservationDir:
            conservationDir = os.path.join(self.outdir, "..", "conservation")
            writeTrackDb_conservation(f, currgenome, conservationDir)

        if self.options.rmskdir:
            writeTrackDb_rmsk(f, os.path.join(self.options.rmskdir, currgenome), self.outdir)

        genomes = self.genomes
        treeGenomes = self.genomes
        
        #Get the neighboring genomes:
        subgenomes = getNeighbors(self.options.tree, currgenome)
        
        #Non composite bed tracks:
        for bigbeddir in self.options.bigbeddirs2:
            writeTrackDb_bigbeds(f, bigbeddir, genomes, currgenome, self.options.properName, False, self.options.tabbed)

        #Composite tracks:
        writeTrackDb_compositeStart(f, self.options.shortLabel, self.options.longLabel, self.options.bigbeddirs, self.options.bigwigdirs, treeGenomes, self.options.properName, self.options.url, self.options.treeFig)
        #if self.options.treeFig:
        #    writeTrackDb_composite_html(os.path.join(self.outdir, "hubCentral.html"), self.options.treeFig)
        for bigbeddir in self.options.bigbeddirs:
            if hasFiles(currgenome, bigbeddir, "bb"):
                writeTrackDb_compositeSubTrack(f, os.path.basename(bigbeddir.rstrip("/")), "dense")
                writeTrackDb_bigbeds(f, bigbeddir, genomes, currgenome, self.options.properName, True, self.options.tabbed)

        for bigwigdir in self.options.bigwigdirs:
            if hasFiles(currgenome, bigwigdir, "bw"):
                writeTrackDb_compositeSubTrack(f, os.path.basename(bigwigdir.rstrip("/")), "dense")
                writeTrackDb_bigwigs(f, bigwigdir, genomes, currgenome, self.options.properName)
                #writeTrackDb_bigwigs(f, bigwigdir, genomes, subgenomes, currgenome, self.options.properName)

        writeTrackDb_compositeSubTrack(f, "Alignments", "full")
        writeTrackDb_snakes(f, self.halfile, genomes, subgenomes, currgenome, self.options.properName, self.options.snpwidth, self.options.selfAlignmentTrack)
        f.close()

############################ UTILITIES FUNCTIONS ###################
def hasFiles(genome, indir, ext):
    for d in os.listdir(indir):
        filepath = os.path.join(indir, d, "%s.%s" %(genome, ext))
        if os.path.exists(filepath):
            return True
    return False

def getLongestSeq(seq2len):
    seqs = sorted( [(seq, len) for seq, len in seq2len.items()], key=lambda item:item[1], reverse=True )
    return seqs[0]

def getGenomeSequencesFromHal(halfile, genome, ucscNames):
    statsfile = getTempFile("%s-seqStats.txt" %genome)
    system("halStats --sequenceStats %s %s > %s" %(genome, halfile, statsfile))

    seq2len = {}
    f = open(statsfile, 'r')
    for line in f:
        if len(line) < 2 or re.search("SequenceName", line):
            continue
        items = line.strip().split(", ")
        seq = items[0]
        if ucscNames:
            seq = seq.split('.')[-1]
        l = int(items[1])
        seq2len[seq] = l
    f.close()
    system("rm %s" %statsfile)

    return seq2len

def getGenomeSequences(halfile, genomes, ucscNames):
    genome2seq2len = {}
    for genome in genomes:
        seq2len = getGenomeSequencesFromHal(halfile, genome, ucscNames)
        if len(seq2len) == 0:
            sys.stderr.write("Warning: genome %s contains 0 sequence - no browser was made.\n" %genome)
        else:
            genome2seq2len[genome] = seq2len
    return genome2seq2len

def getChromSizesFromHal(halfile, genome, outfile):
    system("halStats --chromSizes %s %s > %s" % (genome, halfile, outfile))

def getChromSizes(halfile, seq2len, outfile):
    f = open(outfile, 'w')
    for s, l in seq2len.items():
        if l > 0:
            f.write("%s\t%d\n" %(s, l))
    f.close()

def getGenomesFromHal(halfile):
    #Get a list of all genomes from the output of halStats
    statsfile = getTempFile("halStats.txt")
    system("halStats --genomes %s > %s" %(halfile, statsfile))
    
    f = open(statsfile, 'r')
    genomes = f.readline().strip().split()
    f.close()

    #clean up
    system("rm %s" %statsfile)
    
    return genomes

def linkTwoBitSeqFile(genome, twobitdir, outdir):
    twobitfile = os.path.join(outdir, "%s.2bit" %genome)
    intwobitfile = os.path.abspath( os.path.join(twobitdir, "%s.2bit" %genome) )
    if not os.path.exists(twobitfile):
        system("ln -s %s %s" %(intwobitfile, twobitfile))

def addOptions(parser):
    parser.add_argument('--cpHalFileToOut', dest='cpHal', action='store_true', default=False, help='If specified, copy the input halfile to the output directory (instead of just make a softlink). Default=False')
    parser.add_argument('--ucscNames', dest='ucscNames', action='store_true',
                      default=False,
                      help='Assume that sequence headers use the UCSC '
                      'naming convention, (i.e. "genome.chr"), and  '
                      'attempt to rename the sequences so that their names '
                      'will end up as "chr"')
    addHubOptions(parser)
    addLodOptions(parser)
    addBedOptions(parser)
    addWigOptions(parser)
    addRmskOptions(parser)
    addGcOptions(parser)
    addAlignabilityOptions(parser)
    addConservationOptions(parser)
    addExclusiveRegionOptions(parser)
    addSnakeOptions(parser)

def checkOptions(parser, options):

    # This is not a full Toil port.  Files will still be accessed
    # directly from disk
    options.halfile = os.path.abspath(options.halfile)
    options.outputDirectory = os.path.abspath(options.outputDirectory)

    jobStoreType, locator = Toil.parseLocator(options.jobStore)
    if jobStoreType != "file":
        raise RuntimeError("only local jobStores are supported")
    
    if not os.path.exists(options.halfile):
        raise RuntimeError("Input hal file %s does not exist.\n" %options.halfile)
    if not os.path.exists(options.outputDirectory):
        system("mkdir -p %s" %options.outputDirectory)
    elif not os.path.isdir(options.outputDirectory):
        raise RuntimeError("Output directory specified (%s) is not a directory\n" %options.outputDirectory)
    
    options.snpwidth = None
    checkHubOptions(parser, options)
    checkBedOptions(parser, options)
    checkWigOptions(parser, options)
    checkRmskOptions(parser, options)
    checkConservationOptions(parser, options)

def main():
    usage = '%prog <halFile> <outputDirectory> [options]'
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('halfile', help='input hal file')
    parser.add_argument('outputDirectory', help='output directory')
    addOptions(parser)
    options = parser.parse_args()
    checkOptions(parser, options)

    with Toil(options) as toil:
        toil.start(Setup(options.halfile, options.outputDirectory, options))
        
if __name__ == '__main__':
    main()

