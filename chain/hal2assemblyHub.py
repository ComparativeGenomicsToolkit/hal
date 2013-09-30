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
from sonLib.bioio import system  #(may need to remove this dependency on sonLib)
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

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
                self.addChildTarget( GetCGpercent(genomedir, genome) ) #genomedir/genome.cg.bw
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
                self.addChildTarget( LiftoverBedFiles(beddir, self.halfile, self.genome2seq2len, bigbeddir, self.outdir) )
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

class GetCGpercent( Target ):
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

class GetAlignability( Target ):
    def __init__(self, genomedir, genome, halfile):
        Target.__init__(self)
        self.genomedir = genomedir
        self.genome = genome
        self.halfile = halfile

    def run(self):
        outfile = os.path.join(self.genomedir, "%s.alignability.bw" %self.genome)
        tempwig = os.path.join(self.genomedir, "%s.alignability.wig" %self.genome)
        system("halAlignability %s %s > %s" %(self.halfile, self.genome, tempwig))
        chrsizefile = os.path.join(self.genomedir, "chrom.sizes")
        system("wigToBigWig %s %s %s" %(tempwig, chrsizefile, outfile))
        system("rm -f %s" %tempwig)

class GetConservationFiles( Target ):
    '''Compute conservation bigWig for the input genome and lift it over to all other genomes
    '''
    def __init__(self, halfile, outdir, options):
        Target.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self):
        #First create the model:
        #halPhyloPTrain.py ../../out.hal reference genesMin66Samples.bed ref.mod --numProc 16 --tree rootedTree.nw
        modfile = os.path.join(self.outdir, "%s.mod" % self.options.conservationGenomeName)
        cmd = "halPhyloPTrain.py %s %s %s %s --numProc %d" %(self.halfile, \
                self.options.conservationGenomeName, self.options.conservation, \
                modfile, self.options.conservationNumProc)
        if self.options.conservationTree:
            cmd += " --tree %s" %self.options.conservationTree
        system(cmd)

        self.setFollowOnTarget( GetConservationFiles2(self.halfile, self.outdir, modfile, self.options.conservationNumProc) )

class GetConservationFiles2( Target ):
    '''Modify mod file to convert small branch lengths (change all the xxxe-1y to xxxe-08)
       Then do liftover
    '''
    def __init__(self, halfile, outdir, modfile, numproc):
        Target.__init__(self)
        self.halfile = halfile
        self.outdir = outdir
        self.modfile = modfile
        self.numproc = numproc

    def run(self):
        newmodfile = "%s-modified" %self.modfile
        #modify small branch lengths (change all the xxxe-1y to xxxe-10)
        system("sed 's/e-1./e-08/g' %s > %s" %(self.modfile, newmodfile))
        #get conservation bigwig and liftover files:
        cmd = "halTreePhyloP.py %s %s %s --bigWig --numProc %d" %(self.halfile, newmodfile, self.outdir, self.numproc)
        system(cmd)

class LiftoverBedFiles( Target ):
    def __init__(self, indir, halfile, genome2seq2len, bigbeddir, outdir):
        Target.__init__(self)
        self.indir = indir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.bigbeddir = bigbeddir
        self.outdir = outdir

    def run(self):
        #beddir has the hierachy: indir/genome/chr1.bed, chr2.bed...
        #for each genome in beddir, lifeover the bed records of that genome to the coordinate of all other genomes
         
        #liftover bed file of each genome with available beds to all genomes
        genomes = self.genome2seq2len.keys()
        tempbeds = []
        
        for genome in os.listdir(self.indir):
            if genome not in genomes:
                continue
            genomeindir = os.path.join(self.indir, genome)
            assert os.path.isdir(genomeindir)

            #Create bed directory for current genome
            genomeoutdir = os.path.join(self.bigbeddir, genome)
            system("mkdir -p %s" %genomeoutdir)
        
            #get all the bed files (".bed" ext) and as files if available (".as" ext) 
            bedfiles, asfile, extrafields, numfield = readBedDir(genomeindir)

            #Copy as file to bigbed dir:
            if asfile:
                system("cp %s %s" %(asfile, os.path.join(genomeoutdir, "%s.as" %genome)))
            elif numfield > 12: #does not have .as file, and have more than 12 fields, just treat as 12 fields
                numfield = 12

            #Concatenate all the input bed files and convert it into bigbed to outdir/genome/genome.bb
            tempbed = "%s-temp.bed" % os.path.join(genomeoutdir, genome)
            system( "cat %s/*bed | cut -f-%d > %s" %(genomeindir, numfield, tempbed) )
            system( "bedSort %s %s" % (tempbed, tempbed) )

            outbigbed = os.path.join(genomeoutdir, "%s.bb" %genome) 
            chrsizefile = os.path.join(self.outdir, genome, "chrom.sizes")
            if not asfile:
                system( "bedToBigBed -type=bed%d -extraIndex=name %s %s %s" %(numfield, tempbed, chrsizefile, outbigbed) )
            else:
                numextra = len(extrafields)
                if numextra > 0:
                    type="bed%d+%d" %(numfield - numextra, numextra)
                    extraIndex = "name,%s" % ",".join(extrafields)
                else:
                    type="bed%d" %numfield
                    extraIndex = "name"
                system( "bedToBigBed -as=%s -type=%s -extraIndex=%s %s %s %s" %(asfile, type, extraIndex, tempbed, chrsizefile, outbigbed) )

            #Liftover to all other genomes:
            for othergenome in genomes:
                if othergenome == genome:
                    continue
                self.addChildTarget( LiftoverBed(genomeoutdir, tempbed, asfile, extrafields, numfield, genome, othergenome, self.halfile, self.outdir) )
            tempbeds.append( tempbed )
        self.setFollowOnTarget( CleanupLiftoverBedFiles(tempbeds) )

class CleanupLiftoverBedFiles(Target):
    def __init__(self, files):
        Target.__init__(self)
        self.files = files

    def run(self):
        if len(self.files) > 0:
            system( "rm %s" % " ".join(self.files) ) #cleanup

class LiftoverBed( Target ):
    def __init__(self, genomeoutdir, bed, asfile, extrafields, numfield, genome, othergenome, halfile, outdir):
        Target.__init__(self)
        self.genomeoutdir = genomeoutdir
        self.bed = bed
        self.asfile = asfile
        self.extrafields = extrafields
        self.numfield = numfield
        self.genome = genome
        self.othergenome = othergenome
        self.halfile = halfile
        self.outdir = outdir

    def run(self):
        liftovertempbed = "%s.bed" % os.path.join(self.genomeoutdir, self.othergenome)
        if len(self.extrafields) > 0:
            system("halLiftover %s %s %s %s %s --keepExtra" %(self.halfile, self.genome, self.bed, self.othergenome, liftovertempbed))
        else:
            system("halLiftover %s %s %s %s %s" %(self.halfile, self.genome, self.bed, self.othergenome, liftovertempbed))
        
        system("bedSort %s %s" %(liftovertempbed, liftovertempbed))
        outbigbed = os.path.join(self.genomeoutdir, "%s.bb" %self.othergenome)
        chrsizefile = os.path.join(self.outdir, self.othergenome, "chrom.sizes")
        if os.stat(liftovertempbed).st_size > 0:#make sure the file is not empty
            if not self.asfile:
                #system("bedToBigBed %s %s %s" %(liftovertempbed, chrsizefile, outbigbed))
                system("bedToBigBed -extraIndex=name %s %s %s" %(liftovertempbed, chrsizefile, outbigbed))
                ##system( "bedToBigBed -type=bed%d -extraIndex=name %s %s %s" %(numfield, tempbed, chrsizefile, outbigbed) )
            else:
                numextra = len(self.extrafields)
                if numextra > 0:
                    type="bed%d+%d" %(self.numfield - numextra, numextra)
                    extraIndex = "name,%s" % ",".join(self.extrafields)
                else:
                    type="bed%d" %self.numfield
                    extraIndex = "name"
                #system( "bedToBigBed -as=%s -type=%s %s %s %s" %(self.asfile, type, liftovertempbed, chrsizefile, outbigbed) )
                system( "bedToBigBed -as=%s -type=%s -extraIndex=%s %s %s %s" %(self.asfile, type, extraIndex, liftovertempbed, chrsizefile, outbigbed) )

        #Cleanup:
        system("rm %s" % liftovertempbed)

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
        lodtxtfile = ''
        loddir = ''
        options.lodOpts = ''
        if options.lodMaxBlock is not None:
            options.lodOpts += '--maxBlock %d ' % options.lodMaxBlock
        if options.lodScale is not None:
            options.lodOpts += '--scale %f ' % options.lodScale
        if options.lodMaxDNA is not None:
            options.lodOpts += '--maxDNA %d ' % options.lodMaxDNA
        if options.lodInMemory is True:
            options.lodOpts += '--inMemory '
        if options.lodNumProc is not None:
            options.lodOpts += '--numProc %d ' % options.lodNumProc
        if options.lodMinSeqFrac is not None:
            options.lodOpts += '--minSeqFrac %f ' % options.lodMinSeqFrac
        if options.lodChunk is not None:
            options.lodOpts += '--chunk %d ' % options.lodChunk
        if len(options.lodOpts) > 0:
            options.lod = True
        if options.lod:
            lodtxtfile, loddir = getLodFiles(localHalfile, options, self.outdir)
        
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
            writeTrackDb_cgPercent(f, currgenome)
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

def writeTrackDb_cgPercent(f, genome):
    f.write("track cgPercent\n")
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
    f.write("\n")

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
    f.write("autoScale Off\n")
    f.write("maxHeightPixels 128:36:16\n")
    f.write("graphTypeDefault Bar\n")
    f.write("gridDefault OFF\n")
    f.write("color 0,0,0\n")
    f.write("altColor 128,128,128\n")
    f.write("viewLimits 0:%d\n" %genomeCount)
    f.write("\n")

def writeTrackDb_conservation(f, genome, conservationDir):
    wigfile = os.path.join(conservationDir, "%s_phyloP.bw" %genome)
    if os.path.exists(wigfile):
        f.write("track conservation\n")
        f.write("longLabel Conservation\n")
        f.write("shortLabel Conservation\n")
        f.write("type bigWig -1 1\n")
        f.write("group map\n")
        f.write("visibility dense\n")
        f.write("windowingFunction Mean\n")
        f.write("bigDataUrl ../conservation/%s_phyloP.bw\n" %genome)
        
        f.write("priority 2\n")
        f.write("autoScale Off\n")
        f.write("maxHeightPixels 128:36:16\n")
        f.write("graphTypeDefault Bar\n")
        f.write("gridDefault OFF\n")
        f.write("color 0,0,0\n")
        f.write("altColor 128,128,128\n")
        f.write("viewLimits -1:1\n")
        f.write("\n")

def writeTrackDb_rmsk(f, rmskdir, genomedir):
    if not os.path.exists(rmskdir):
        return
    f.write("track repeatMasker_\n")
    f.write("compositeTrack on\n")
    f.write("shortLabel RepeatMasker\n")
    f.write("longLabel Repeating Elements by RepeatMasker\n")
    f.write("group map\n")
    f.write("visibility dense\n")
    f.write("type bed 3 .\n")
    f.write("noInherit on\n")
    f.write("\n")
    
    system("ln -s %s %s" %(os.path.abspath(rmskdir), os.path.join(genomedir, "repeatMasker")))
    files = os.listdir(rmskdir)
    for i, file in enumerate(files):
        element = file.split('.')[0]
        f.write("\ttrack repeatMasker%s\n" %element)
        f.write("\tparent repeatMasker_\n")
        f.write("\tshortLabel %s\n" %element)
        f.write("\tlongLabel %s Repeating Elements by RepeatMasker\n" %element)
        f.write("\tpriority %d\n" %i)
        f.write("\tspectrum on\n")
        f.write("\tmaxWindowToDraw 10000000\n")
        f.write("\tcolorByStrand 50,50,150 150,50,50\n")
        f.write("\ttype bigBed 6 +\n")
        f.write("\tbigDataUrl repeatMasker/%s.bb\n" %element)
        f.write("\n")
    f.write("\n")

def writeTrackDb_bigbeds(f, bigbeddir, genomes, currgenome, properName):
    annotation = os.path.basename(bigbeddir)
    genome2priority = {}
    for i, genome in enumerate(genomes):
        if genome == currgenome:
            genome2priority[genome] = 1
        else:
            genome2priority[genome] = i + 2
    
    for genome in os.listdir(bigbeddir):
        bbfile = os.path.join(bigbeddir, genome, "%s.bb" %currgenome)
        if not os.path.exists(bbfile):
            continue
        #see if there is an as file:
        searchIndexStr = "name"
        asfile = os.path.join(bigbeddir, genome, "%s.as" %genome)
        if os.path.exists(asfile):
            extrafields, numfield = getBedExtraFieldsFromAsFile(asfile)
            fields = "%d" %numfield
            if len(extrafields) > 0:
                fields = "%d +" %(numfield - len(extrafields))
                searchIndexStr += ",%s" %(",".join(extrafields))
            else:
                fields += " ."
        else:
            numfield = getBedNumField(bbfile)
            fields = "%d ." %numfield
        #start writing track
        genomeProperName = genome
        if genome in properName:
            genomeProperName = properName[genome]
        priority = 1
        if genome in genome2priority:
            priority = genome2priority[genome]

        f.write("track %s%s\n" % (annotation, genome))
        if genome == currgenome:
            f.write("longLabel %s %s\n" % (genomeProperName, annotation))
        else:
            f.write("longLabel %s Lifted-over %s\n" % (genomeProperName, annotation))
        f.write("priority %d\n" %priority)
        f.write("shortLabel %s%s\n" % (genomeProperName, annotation))
        f.write("bigDataUrl ../liftoverbeds/%s\n" % os.path.join( annotation, genome, "%s.bb" % currgenome ) )
        f.write("type bigBed %s\n" %fields)
        f.write("group annotation%s\n" %annotation)
        if numfield >=4:
        #if numfield >=4 and genome == currgenome: #DEBUG
            f.write("searchIndex %s\n" %searchIndexStr)
        #if not re.search('gene', annotation): #HACK
        f.write("itemRgb On\n")
        if genome == currgenome or not re.search('Gene', annotation): #HACK
            f.write("visibility dense\n")
        else:
            f.write("visibility hide\n")
        f.write("\n")

def writeTrackDb_snakes(f, halfile, genomes, currgenome, properName):
    for i, genome in enumerate(genomes):
        if re.search(genome, currgenome): #current genome
            continue
        #SNAKE TRACKS
        genomeProperName = genome
        if genome in properName:
            genomeProperName = properName[genome]
        f.write("track snake%s\n" %genome)
        f.write("longLabel %s\n" %genomeProperName)
        f.write("shortLabel %s\n" %genomeProperName)
        f.write("otherSpecies %s\n" %genome)
        f.write("visibility full\n")
        f.write("priority %d\n" %(i + 2))
        f.write("bigDataUrl %s\n" % halfile)
        f.write("type halSnake\n")
        f.write("group snake\n")
        f.write("\n")

def writeDescriptionFile(genome, outdir):
    filename = os.path.join(outdir, "description.html")
    f = open(filename, 'w')
    f.write("%s\n" %genome)
    f.close()
    return

def fixLodFilePath(lodtxtfile, localHalfile, outdir):
    #fix the path of the original hal file to point to the created
    #link relative to the output directory
    relPath = os.path.relpath(localHalfile, start=outdir)
    lodTxtBuf = ''
    for line in open(lodtxtfile):
        tokens = line.split()
        if len(tokens) == 2 and tokens[0] == '0':
            lodTxtBuf += '0 %s\n' % relPath
        else:
            lodTxtBuf += line
    with open(lodtxtfile, 'w') as lodFile:
        lodFile.write(lodTxtBuf)
    
def getLodFiles(localHalfile, options, outdir):
    lodtxtfile = os.path.join(outdir, "lod.txt") #outdir/lod.txt
    loddir = os.path.join(outdir, "lod") #outdir/lod
    if options.lodtxtfile and options.loddir: #if lod files were given, then just make soft links to them
        if os.path.exists(lodtxtfile):
            if os.path.abspath(lodtxtfile) != os.path.abspath(options.lodtxtfile):
                system("rm %s" %lodtxtfile)
                system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))

        if os.path.exists(loddir):
            if os.path.abspath(loddir) != os.path.abspath(options.loddir):
                if os.path.islink(loddir):
                    system("rm %s" %loddir)
                else:
                    system("rm -Rf %s" %loddir)
                loddir = os.path.join(outdir, os.path.basename(options.loddir))
                system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
    else: #if lod files were not given, create them using halLodInterpolate.py
        system("halLodInterpolate.py %s %s --outHalDir %s %s" %(localHalfile, lodtxtfile, loddir, options.lodOpts))
        fixLodFilePath(lodtxtfile, localHalfile, outdir)
    return lodtxtfile, loddir

def writeGroupFile(outdir, annotations):
    filename = os.path.join(outdir, "groups.txt")
    f = open(filename, 'w')
    f.write("name user\n")
    f.write("label Custom\n")
    f.write("priority 1\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")

    f.write("name map\n")
    f.write("label Mapping\n")
    f.write("priority 2\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")

    f.write("name snake\n")
    f.write("label Alignment Snakes\n")
    f.write("priority 3\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")
    
    for annotation in annotations:
        f.write("name annotation%s\n" %annotation)
        f.write("label %s Annotations\n" % annotation.capitalize() )
        f.write("priority 3\n")
        f.write("defaultIsClosed 1\n")
        f.write("\n")   

    f.write("name exp\n")
    f.write("label Experimental\n")
    f.write("priority 4\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")
    
    f.close()

def writeHubFile(outdir, options):
    hubfile = os.path.join(outdir, "hub.txt")
    f = open(hubfile, "w")
    f.write("hub %s\n" %options.hubLabel)
    f.write("shortLabel %s\n" %options.shortLabel)
    f.write("longLabel %s\n" %options.longLabel)
    f.write("genomesFile genomes.txt\n")
    f.write("email %s\n" %options.email)
    f.close()

def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        items.append(line.strip())
    f.close()
    return items

def readRename(file):
    name2new = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        if len(items) >=2:
            name2new[items[0]] = items[1]
    f.close()
    return name2new

def getFilesByExt(indir, fileExtension):
    #Return all files in indir that have "fileExtension"
    allfiles = os.listdir(indir)
    files = []
    for f in allfiles:
        if os.path.splitext(f) == ".%s" %fileExtension :
            files.append(f)
    return files

def readBedDir(indir):
    #Get all the bed files in "indir", as well as the accompanying ".as" files if present
    bedfiles = []
    asfile = None 
    extrafields = []
    numfield = 0
    allfiles = os.listdir(indir)

    for f in allfiles:
        ext = os.path.splitext(f)[-1]
        if ext == ".bed":
            if os.stat(os.path.join(indir, f)).st_size > 0:#make sure the file is not empty
                bedfiles.append(f)
        elif ext == ".as":
            asfile = os.path.join(indir, f)

    #Check to make sure asfile and bedfiles have the same number of fields:
    if asfile:
        extrafields, numfield = getBedExtraFieldsFromAsFile(asfile)
    elif len(bedfiles) > 0:
        numfield = getFileColumnCount(os.path.join(indir, bedfiles[0]))
    
    for b in bedfiles:
        bedfile = os.path.join(indir, b)
        assert getFileColumnCount(bedfile) == numfield

    return bedfiles, asfile, extrafields, numfield

def getFileColumnCount(file):
    f = open(file, 'r')
    firstline = f.readline()
    items = firstline.split()
    f.close()
    return len(items)

def getBedNumField(bbfile): 
    tempbed = "%s-TEMP-%s" %(bbfile, time.time())
    system("bigBedToBed %s %s" %(bbfile, tempbed))
    numfield = getFileColumnCount(tempbed)
    system("rm %s" %tempbed)
    return numfield

def getBedExtraFieldsFromAsFile(asfile):
    numfield = 0
    fields = []
    standardFields = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts'] #standard fields of a bed entry
    
    f = open(asfile, 'r')
    begin = False
    for line in f:
        line = line.strip()
        if re.search("^\(", line):
            begin = True
            line = line.lstrip("(")
        if len(line) == 0 or not begin or re.search("\)", line):
            continue
        items = line.split('\t')
        assert len(items) >= 2
        field = items[1].rstrip(';')
        numfield += 1
        if field not in standardFields:
            fields.append(field)
    f.close()
    return fields, numfield

def addOptions(parser):
    parser.add_option('--lod', dest='lod', action="store_true", default=False, help='If specified, create "level of detail" (lod) hal files and will put the lod.txt at the bigUrl instead of the original hal file. Default=%default')
    parser.add_option('--lodTxtFile', dest='lodtxtfile', help='"hal Level of detail" lod text file. If specified, will put this at the bigUrl instead of the hal file. Default=%default')
    parser.add_option('--lodDir', dest='loddir', help='"hal Level of detail" lod dir. If specified, will put this at the bigUrl instead of the hal file. Default=%default')
    parser.add_option('--lodMaxBlock', dest='lodMaxBlock', type='int', help='Maximum number of blocks to display in a hal level of detail. Default=%default', default=None)
    parser.add_option('--lodScale', dest='lodScale', type='float', help='Scaling factor between two successive levels of detail. Default=%default.', default=None)
    parser.add_option('--lodMaxDNA', dest='lodMaxDNA', type='int', help='Maximum query length that will such that its hal level of detail will contain nucleotide information. Default=%default.', default=None)
    parser.add_option('--lodInMemory', dest='lodInMemory', action='store_true', help='Load entire hal file into memory when generating levels of detail instead of using hdf5 cache. Default=%default.', default=False)
    parser.add_option('--lodNumProc', dest='lodNumProc', type='int', help='Number of levels of detail to generate concurrently in parallel processes', default=None)
    parser.add_option('--lodMinSeqFrac', dest='lodMinSeqFrac', type='float', help='Minumum sequence length to sample as fraction of step size for level of detail generation: ie sequences with length <= floor(minSeqFrac * step) are ignored. Use default from halLodExtract if not set.', default=None)
    parser.add_option('--lodChunk', dest='lodChunk', type='int', help='HDF5 chunk size for generated levels of detail.', default=None)
    parser.add_option('--cpHalFileToOut', dest='cpHal', action='store_true', default=False, help='If specified, copy the input halfile to the output directory (instead of just make a softlink). Default=%default')
    parser.add_option('--hub', dest='hubLabel', default='myHub', help='a single-word name of the directory containing the track hub files. Not displayed to hub users. Default=%default')
    parser.add_option('--shortLabel', dest='shortLabel', default='my hub', help='the short name for the track hub. Suggested maximum length is 17 characters. Displayed as the hub name on the Track Hubs page and the track group name on the browser tracks page. Default=%default')
    parser.add_option('--longLabel', dest='longLabel', default='my hub', help='a longer descriptive label for the track hub. Suggested maximum length is 80 characters. Displayed in the description field on the Track Hubs page. Default=%default')
    parser.add_option('--email', dest='email', default='NoEmail', help='the contact to whom questions regarding the track hub should be directed. Default=%default')
    parser.add_option('--bedDirs', dest='beddirs', help='comma separated list of directories containing bed files of the input genomes. Each directory represents a type of annotation. The annotations of each genome will then be liftovered to all other genomes in the MSA. Example: "genes,genomicIsland,tRNA". Format of each directory: bedDir/ then genome1/ then chr1.bed, chr2.bed... Default=%default' )
    parser.add_option('--finalBigBedDirs', dest='bbdirs', help='comma separated list of directories containing final big bed files to be displayed. No liftover will be done for these files. Each directory represents a type of annotation. Example: "genes,genomicIsland,tRNA". Format of each directory: bbDir/ then queryGenome/ then targetGenome1.bb, targetGenome2.bb ... (so annotation of queryGenome has been mapped to targetGenomes and will be display on the targetGenome browsers). Default=%default' )
    parser.add_option('--rmskDir', dest='rmskdir', help="Directory containing repeatMasker's output files for each genome. Format: rmskDir/ then genome1/ then genome.rmsk.SINE.bb, genome.rmsk.LINE.bb, ... Default=%default")
    parser.add_option('--gcContent', dest='gcContent', action='store_true', default=False, help='If specified, make GC-content tracks. Default=%default')
    parser.add_option('--alignability', dest='alignability', action='store_true', default=False, help='If specified, make Alignability tracks. Default=%default')
    parser.add_option('--genomes', dest='genomes', help='File specified list of genomes to make browser for. If specified, only create browsers for these genomes in the order provided by the list. Otherwise create browsers for all genomes in the input hal file')
    parser.add_option('--conservation', dest='conservation', help='Bed file providing regions to calculate the conservation tracks.')
    parser.add_option('--conservationGenomeName', dest='conservationGenomeName', help='Name of the genome of the bed file provided in the "--conversation" option')
    parser.add_option('--conservationTree', dest='conservationTree', help='Optional. Newick tree for the conservation track')
    parser.add_option('--conservationNumProc', dest='conservationNumProc', type='int', default=1, help='Optional. Number of processors to run conservation')
    parser.add_option('--conservationDir', dest='conservationDir', help='Directory contains conservation bigwigs')
    parser.add_option('--rename', dest='rename', help='File that maps halfile genomeNames to names displayed on the browser. Format: <halGenomeName>\\t<genomeNameToDisplayOnBrowser>. Default=%default') 

def checkOptions(parser, args, options):
    if len(args) < 2:
        parser.error("Required two input arguments, %d was provided\n" %len(args))
    if not os.path.exists(args[0]):
        parser.error("Input hal file %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        system("mkdir -p %s" %args[1])
    elif not os.path.isdir(args[1]):
        parser.error("Output directory specified (%s) is not a directory\n" %args[1])
    if options.genomes:
        options.genomes = readList(options.genomes)
    options.properName = {}
    if options.rename and os.path.exists(options.rename):
        options.properName = readRename(options.rename)

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
    from hal.chain.hal2assemblyHub import *
    main()

