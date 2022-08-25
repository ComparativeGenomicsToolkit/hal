#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Creating bed tracks and lifted-over bed tracks for the hubs
"""
import os, re, time
from sonLib.bioio import system  
from toil.job import Job
from optparse import OptionGroup
from hal.assemblyHub.assemblyHubCommon import CleanupFiles, getProperName
from hal.assemblyHub.bedCommon import filterLongIntrons, tabifyBed, untabifyBed

class LiftoverBedFiles( Job ):
    def __init__(self, indir, halfile, genome2seq2len, bigbeddir, noLiftover, tab, outdir, options):
        Job.__init__(self)
        self.indir = indir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.bigbeddir = bigbeddir
        self.noLiftover = noLiftover
        self.tab = tab
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        #beddir has the hierachy: indir/genome/chr1.bed, chr2.bed...
        #for each genome in beddir, lifeover the bed records of that genome to the coordinate of all other genomes
         
        #liftover bed file of each genome with available beds to all genomes
        genomes = list(self.genome2seq2len.keys())
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
            bedfiles, asfile, extrafields, numfield = readBedDir(genomeindir, self.tab)
            if numfield < 3:
                # This is an empty (probably from an automated
                # process) or otherwise malformed bed. Whine to the
                # user and then attempt to go as far as possible
                # anyway.
                self.logToMaster("WARNING: input bed files in %s have less "
                                 "than 3 fields, or are completely empty. "
                                 "Proceeding anyway." % genomeindir)
                numfield = 3

            #Copy as file to bigbed dir:
            if asfile:
                system("cp %s %s" %(asfile, os.path.join(genomeoutdir, "%s.as" %genome)))
            elif numfield > 12: #does not have .as file, and have more than 12 fields, just treat as 12 fields
                numfield = 12

            #Concatenate all the input bed files and convert it into bigbed to outdir/genome/genome.bb
            tempbed = "%s-temp.bed" % os.path.join(genomeoutdir, genome)
            system( "cat %s/*bed | cut -f-%d > %s" %(genomeindir, numfield, tempbed) )
            #system( "bedSort %s %s" % (tempbed, tempbed) )
            filterbed = "%s-temp-filtered.bed" %os.path.join(genomeoutdir, genome)
            filterLongIntrons(tempbed, filterbed, 100000, self.tab,
                              self.options.ucscNames)
            # bedSort expects tab-separated beds, so we have to do some
            # format gymnastics here.
            if not self.tab:
                tabifyBed(filterbed)
            system( "bedSort %s %s" % (filterbed, tempbed) )
            if not self.tab:
                untabifyBed(tempbed)

            outbigbed = os.path.join(genomeoutdir, "%s.bb" %genome) 
            chrsizefile = os.path.join(self.outdir, genome, "chrom.sizes")
            if not asfile:
                # Index on the 'name' field if the bed has one
                indexParameter = "-extraIndex=name" if numfield >= 4 else ""
                cmd = "bedToBigBed -type=bed%d %s %s %s %s" %(numfield, indexParameter, tempbed, chrsizefile, outbigbed)
                if self.tab:
                    cmd = "bedToBigBed -tab -type=bed%d %s %s %s %s" %(numfield, indexParameter, tempbed, chrsizefile, outbigbed)
                system( cmd )
            else:
                assert numfield >= 4 # -extraIndex=name will fail if this is not true.
                numextra = len(extrafields)
                if numextra > 0:
                    type="bed%d+%d" %(numfield - numextra, numextra)
                    extraIndex = "name,%s" % ",".join(extrafields)
                else:
                    type="bed%d" %numfield
                    extraIndex = "name"
                cmd = "bedToBigBed -as=%s -type=%s -extraIndex=%s %s %s %s" %(asfile, type, extraIndex, tempbed, chrsizefile, outbigbed)
                if self.tab:
                    cmd = "bedToBigBed -tab -as=%s -type=%s -extraIndex=%s %s %s %s" %(asfile, type, extraIndex, tempbed, chrsizefile, outbigbed)
                system( cmd )

            #Liftover to all other genomes:
            if not self.noLiftover:
                for othergenome in genomes:
                    if othergenome == genome:
                        continue
                    self.addChild( LiftoverBed(genomeoutdir, tempbed, self.tab, asfile, extrafields, numfield, genome, othergenome, self.halfile, self.outdir, self.options) )
            tempbeds.append( tempbed )
            tempbeds.append( filterbed )
        self.addFollowOn( CleanupFiles(tempbeds) )

class LiftoverBed( Job ):
    def __init__(self, genomeoutdir, bed, tab, asfile, extrafields, numfield, genome, othergenome, halfile, outdir, options):
        Job.__init__(self)
        self.genomeoutdir = genomeoutdir
        self.bed = bed
        self.tab = tab
        self.asfile = asfile
        self.extrafields = extrafields
        self.numfield = numfield
        self.genome = genome
        self.othergenome = othergenome
        self.halfile = halfile
        self.outdir = outdir
        self.options = options

    def run(self, fileStore):
        liftovertempbed = "%s.bed" % os.path.join(self.genomeoutdir, self.othergenome)

        # bedSort (and, given that it's --tab option is gone, presumably halLiftover) expects tab-separated beds, so we have to do some
        # format gymnastics here.
        if not self.tab:
            tabifyBed(self.bed)
            
        cmd = "halLiftover %s %s %s %s %s" %(self.halfile, self.genome, self.bed, self.othergenome, liftovertempbed)
        if len(self.extrafields) > 0:
            cmd += " --keepExtra"
        else:
            cmd += " --bedType %d" %self.numfield
        system(cmd) 

        filterbed = "%s-filtered.bed" %os.path.join(self.genomeoutdir, self.othergenome)
        filterLongIntrons(liftovertempbed, filterbed, 100000, self.tab, self.options.ucscNames)
        system( "bedSort %s %s" % (filterbed, liftovertempbed) )
        
        if not self.tab:
            untabifyBed(liftovertempbed)

        outbigbed = os.path.join(self.genomeoutdir, "%s.bb" %self.othergenome)
        chrsizefile = os.path.join(self.outdir, self.othergenome, "chrom.sizes")
        if not self.asfile:
            cmd = "bedToBigBed -type=bed%d %s %s %s" %(self.numfield, liftovertempbed, chrsizefile, outbigbed)
            if self.numfield >= 4:
                cmd += " -extraIndex=name"
        else:
            numextra = len(self.extrafields)
            if numextra > 0:
                type="bed%d+%d" %(self.numfield - numextra, numextra)
                extraIndex = "name,%s" % ",".join(self.extrafields)
            else:
                type="bed%d" %self.numfield
                extraIndex = "name"
            cmd = "bedToBigBed -as=%s -type=%s -extraIndex=%s %s %s %s" %(self.asfile, type, extraIndex, liftovertempbed, chrsizefile, outbigbed)
        if self.tab:
            cmd += " -tab"
        system(cmd)

        #Cleanup:
        system("rm %s" % liftovertempbed)
        system("rm -f %s" % filterbed)

def getPriorities(genomes, currgenome):
    genome2priority = {}
    for i, genome in enumerate(genomes):
        if genome == currgenome:
            genome2priority[genome] = 1
        else:
            genome2priority[genome] = i + 2
    return genome2priority 

def getSearchIndexInfo(bigbeddir, genome, bbfile, tab):
    #see if there is an as file:
    asfile = os.path.join(bigbeddir, genome, "%s.as" %genome)
    searchIndexStr = "name"
    if os.path.exists(asfile):
        extrafields, numfield = getBedExtraFieldsFromAsFile(asfile)
        fields = "%d" %numfield
        if len(extrafields) > 0:
            fields = "%d +" %(numfield - len(extrafields))
            searchIndexStr += ",%s" %(",".join(extrafields))
        else:
            fields += " ."
    else:
        numfield = getBedNumField(bbfile, tab)
        fields = "%d ." %numfield
    return searchIndexStr, fields, numfield

def writeTrackDb_bigbeds(f, bigbeddir, genomes, currgenome, properName, composite, tab):
    annotation = os.path.basename(bigbeddir)
    genome2priority = getPriorities(genomes, currgenome)
    
    for genome in os.listdir(bigbeddir):
        bbfile = os.path.join(bigbeddir, genome, "%s.bb" %currgenome)
        if not os.path.exists(bbfile):
            continue
        searchIndexStr, fields, numfield = getSearchIndexInfo(bigbeddir, genome, bbfile, tab)
        
        #start writing track
        genomeProperName = getProperName(genome, properName)
        priority = 1
        if genome in genome2priority:
            priority = genome2priority[genome]
        indent = ''
        if composite:
            indent = "\t\t"
        f.write("%strack %s%s\n" % (indent, annotation, genome))
        if genome == currgenome:
            f.write("%slongLabel %s %s\n" % (indent, genomeProperName, annotation))
        else:
            f.write("%slongLabel %s Lifted-over %s\n" % (indent, genomeProperName, annotation))
        f.write("%spriority %d\n" %(indent, priority))
        f.write("%stype bigBed %s\n" %(indent, fields))
        f.write("%sgroup annotation%s\n" %(indent, annotation))
        if numfield >=4:
            f.write("%ssearchIndex %s\n" %(indent, searchIndexStr))
        f.write("%sitemRgb On\n" %indent)
        if genome == currgenome:
            f.write("%svisibility dense\n" %indent)
        else:
            f.write("%svisibility hide\n" %indent)
        if composite: #Composite
            f.write("%sshortLabel %s%s\n" % (indent, genomeProperName, annotation))
            f.write("%sbigDataUrl ../liftoverbed/%s\n" % (indent, os.path.join( annotation, genome, "%s.bb" % currgenome )) )
            if genome == currgenome:
                f.write("%sparent hubCentral%s\n" %(indent, annotation))
            else:
                f.write("%sparent hubCentral%s off\n" %(indent, annotation))
            f.write("%ssubGroups view=%s orgs=%s\n" %(indent, annotation, genome))
        else:
            f.write("%sshortLabel %s\n" % (indent, genomeProperName))
            f.write("%sbigDataUrl ../liftoverbed2/%s\n" % (indent, os.path.join( annotation, genome, "%s.bb" % currgenome )) )
        f.write("\n")

def writeTrackDb_bigbeds_hackFakeRow(f, bigbeddir, genomes, currgenome, properName, tab):
    annotation = os.path.basename(bigbeddir)
    genome2priority = getPriorities(genomes, currgenome)
    
    bedGenomes = os.listdir(bigbeddir)#HACK
    bedGenomes2 = os.listdir(bigbeddir)#HACK
    if currgenome not in bedGenomes2:#HACK
        bedGenomes2.append(currgenome)#HACK
    for genome in bedGenomes2:#HACK
        bbfile = os.path.join(bigbeddir, genome, "%s.bb" %currgenome)
        if genome == currgenome and currgenome not in bedGenomes:#HACK
            bbfile = os.path.join(bigbeddir, bedGenomes[1], "%s.bb" %currgenome)#HACK
        if not os.path.exists(bbfile):
            continue
        searchIndexStr, fields, numfield = getSearchIndexInfo(bigbeddir, genome, bbfile, tab)
        
        #start writing track
        genomeProperName = getProperName(genome, properName)
        priority = 1
        if genome in genome2priority:
            priority = genome2priority[genome]

        f.write("\t\ttrack %s%s\n" % (annotation, genome))
        if genome == currgenome:
            if currgenome not in bedGenomes:#HACK
                f.write("\t\tlongLabel %s Lifted-over %s\n" % (properName[bedGenomes[1]], annotation))#HACK
            else:#HACK
                f.write("\t\tlongLabel %s %s\n" % (genomeProperName, annotation))#HACK
        else:
            f.write("\t\tlongLabel %s Lifted-over %s\n" % (genomeProperName, annotation))
        f.write("\t\tpriority %d\n" %priority)
        f.write("\t\tshortLabel %s%s\n" % (genomeProperName, annotation))
        if genome == currgenome and genome not in bedGenomes:#HACK
            f.write("\t\tbigDataUrl ../liftoverbed/%s\n" % os.path.join( annotation, bedGenomes[1], "%s.bb" % currgenome ) )#HACK
        else:#HACK
            f.write("\t\tbigDataUrl ../liftoverbed/%s\n" % os.path.join( annotation, genome, "%s.bb" % currgenome ) )#HACK
        f.write("\t\ttype bigBed %s\n" %fields)
        f.write("\t\tgroup annotation%s\n" %annotation)
        if numfield >=4:
            f.write("\t\tsearchIndex %s\n" %searchIndexStr)
        f.write("\t\titemRgb On\n")
        if genome == currgenome:
            if genome in bedGenomes:
                f.write("\t\tvisibility dense\n")
                f.write("\t\tparent hubCentral%s\n"%annotation)
            else:
                f.write("\t\tvisibility hide\n")
                f.write("\t\tparent hubCentral%s off\n"%annotation)
        else:
            f.write("\t\tvisibility hide\n")
            f.write("\t\tparent hubCentral%s off\n"%annotation)
        f.write("\t\tsubGroups view=%s orgs=%s\n" %(annotation, genome))
        f.write("\n")

def readBedDir(indir, tab):
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
        numfield = getFileColumnCount(os.path.join(indir, bedfiles[0]), tab)
    
    for b in bedfiles:
        bedfile = os.path.join(indir, b)
        assert getFileColumnCount(bedfile, tab) == numfield

    return bedfiles, asfile, extrafields, numfield

def getFileColumnCount(file, tab):
    f = open(file, 'r')
    firstline = f.readline()
    if tab:
        items = firstline.split("\t")
    else:
        items = firstline.split()
    f.close()
    return len(items)

def getBedNumField(bbfile, tab): 
    tempbed = "%s-TEMP-%s" %(bbfile, time.time())
    system("bigBedToBed %s %s" %(bbfile, tempbed))
    numfield = getFileColumnCount(tempbed, tab)
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

#========= OPTIONS ==============
def addBedOptions(parser):
    group = parser.add_argument_group("BED-FORMATTED ANNOTATIONS", "All annotations in bed or bigbed formats.")
    group.add_argument('--bedDirs', dest='beddirs', help='comma separated list of directories containing bed files of the input genomes. Each directory represents a type of annotation. The annotations of each genome will then be liftovered to all other genomes in the MSA. Example: "genes,genomicIsland,tRNA". Format of each directory: bedDir/ then genome1/ then chr1.bed, chr2.bed... ' )
    group.add_argument('--finalBigBedDirs', dest='bbdirs', help='comma separated list of directories containing final big bed files to be displayed. No liftover will be done for these files. Each directory represents a type of annotation. Example: "genes,genomicIsland,tRNA". Format of each directory: bbDir/ then queryGenome/ then targetGenome1.bb, targetGenome2.bb ... (so annotation of queryGenome has been mapped to targetGenomes and will be display on the targetGenome browsers). ' )
    group.add_argument('--bedDirs2', dest='beddirs2', help='Similar to --bedDirs, except these tracks will be kept separately and out of the composite track. ')
    group.add_argument('--finalBigBedDirs2', dest='bbdirs2', help='Similar to --finalBigBedDirs, except these tracks will be kept separately and out of the composite track. ')
    group.add_argument('--noBedLiftover', dest='noBedLiftover', action='store_true', default=False, help='If specified, will not lift over the bed annotations. ')
    group.add_argument('--tabBed', dest='tabbed', action='store_true', default=False, help='If specified, treat tab as the delimiter of all the bed files. Default: any white space.')
    group = parser.add_argument_group(group)

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
    options.bigbeddirs = []
    
    if options.beddirs2:
        dirs = [d.rstrip('/') for d in options.beddirs2.split(',')]
        options.beddirs2 = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Bed directory %s does not exist or is not a directory.\n" %d)
    if options.bbdirs2:
        dirs = [d.rstrip('/') for d in options.bbdirs2.split(',')]
        options.bbdirs2 = dirs
        for d in dirs:
            if not os.path.exists(d) or not os.path.isdir(d):
                parser.error("Bigbed directory %s does not exist or is not a directory.\n" %d)
    options.bigbeddirs2 = []


