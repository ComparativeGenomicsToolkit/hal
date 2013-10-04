#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""Creating bed tracks and lifted-over bed tracks for the hubs
"""
import os, re, time
from sonLib.bioio import system  
from jobTree.scriptTree.target import Target
from optparse import OptionGroup
from hal.assemblyHub.assemblyHubCommon import CleanupFiles

class LiftoverBedFiles( Target ):
    def __init__(self, indir, halfile, genome2seq2len, bigbeddir, noLiftover, outdir):
        Target.__init__(self)
        self.indir = indir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.bigbeddir = bigbeddir
        self.noLiftover = noLiftover
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
            if not self.noLiftover:
                for othergenome in genomes:
                    if othergenome == genome:
                        continue
                    self.addChildTarget( LiftoverBed(genomeoutdir, tempbed, asfile, extrafields, numfield, genome, othergenome, self.halfile, self.outdir) )
            tempbeds.append( tempbed )
        self.setFollowOnTarget( CleanupFiles(tempbeds) )

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

        f.write("\t\ttrack %s%s\n" % (annotation, genome))
        if genome == currgenome:
            f.write("\t\tlongLabel %s %s\n" % (genomeProperName, annotation))
        else:
            f.write("\t\tlongLabel %s Lifted-over %s\n" % (genomeProperName, annotation))
        f.write("\t\tpriority %d\n" %priority)
        f.write("\t\tshortLabel %s%s\n" % (genomeProperName, annotation))
        f.write("\t\tbigDataUrl ../liftoverbeds/%s\n" % os.path.join( annotation, genome, "%s.bb" % currgenome ) )
        f.write("\t\ttype bigBed %s\n" %fields)
        f.write("\t\tgroup annotation%s\n" %annotation)
        if numfield >=4:
        #if numfield >=4 and genome == currgenome: #DEBUG
            f.write("\t\tsearchIndex %s\n" %searchIndexStr)
        #if not re.search('gene', annotation): #HACK
        f.write("\t\titemRgb On\n")
        if genome == currgenome or not re.search('Gene', annotation): #HACK
            f.write("\t\tvisibility dense\n")
        else:
            f.write("\t\tvisibility hide\n")
        f.write("\t\tparent hubCentral%s\n"%annotation)
        f.write("\t\tsubGroups view=%s orgs=%s\n" %(annotation, genome))
        f.write("\n")

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

def addBedOptions(parser):
    group = OptionGroup(parser, "BED-FORMATTED ANNOTATIONS", "All annotations in bed or bigbed formats.")
    group.add_option('--bedDirs', dest='beddirs', help='comma separated list of directories containing bed files of the input genomes. Each directory represents a type of annotation. The annotations of each genome will then be liftovered to all other genomes in the MSA. Example: "genes,genomicIsland,tRNA". Format of each directory: bedDir/ then genome1/ then chr1.bed, chr2.bed... Default=%default' )
    group.add_option('--finalBigBedDirs', dest='bbdirs', help='comma separated list of directories containing final big bed files to be displayed. No liftover will be done for these files. Each directory represents a type of annotation. Example: "genes,genomicIsland,tRNA". Format of each directory: bbDir/ then queryGenome/ then targetGenome1.bb, targetGenome2.bb ... (so annotation of queryGenome has been mapped to targetGenomes and will be display on the targetGenome browsers). Default=%default' )
    group.add_option('--noBedLiftover', dest='noBedLiftover', action='store_true', default=False, help='If specified, will not lift over the bed annotations. Default=%default')
    parser.add_option_group(group)

def checkBedOptions(parser, options):
    if options.beddirs:
        dirs = [d.rstrip('/') for d in options.beddirs.split(',')]
        options.beddirs = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Bed directory %s does not exist or is not a directory.\n" %d)
    if options.bbdirs:
        dirs = [d.rstrip('/') for d in options.bbdirs.split(',')]
        options.bbdirs = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Bigbed directory %s does not exist or is not a directory.\n" %d)


