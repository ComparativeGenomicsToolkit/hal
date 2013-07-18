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

import os, sys, re
from optparse import OptionParser
from sonLib.bioio import system  #(may need to remove this dependency on sonLib)

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

def liftoverBedFiles(indir, halfile, genome2seq2len, bigbeddir, outdir):
    #beddir has the hierachy: indir/genome/chr1.bed, chr2.bed...
    #for each genome in beddir, lifeover the bed records of that genome to the coordinate of all other genomes

    #liftover bed file of each genome with availabel beds to all genomes
    genomes = genome2seq2len.keys()
    for genome in os.listdir(indir):
        if genome not in genomes:
            continue
        genomeindir = os.path.join(indir, genome)
        assert os.path.isdir(genomeindir)

        #Create bed directory for current genome
        genomeoutdir = os.path.join(bigbeddir, genome)
        system("mkdir -p %s" %genomeoutdir)

        #Concatenate all the input bed files and convert it into bigbed to outdir/genome/genome.bb
        tempbed = "%s-temp.bed" % os.path.join(genomeoutdir, genome)
        system( "cat %s/*bed > %s" %(genomeindir, tempbed) )
        system( "bedSort %s %s" % (tempbed, tempbed) )

        outbigbed = os.path.join(genomeoutdir, "%s.bb" %genome) 
        chrsizefile = os.path.join(outdir, genome, "chrom.sizes")
        system( "bedToBigBed %s %s %s" %(tempbed, chrsizefile, outbigbed) )

        #Liftover to all other genomes:
        for othergenome in genomes:
            if othergenome == genome:
                continue
            liftovertempbed = "%s.bed" % os.path.join(genomeoutdir, othergenome)
            system("halLiftover %s %s %s %s %s" %(halfile, genome, tempbed, othergenome, liftovertempbed))
            #if re.search("Anc", othergenome): #HACK
            #    system("awk ' $0 !~ /#/ {split($1, arr, \".\"); print arr[2] \"\t\" $2 \"\t\" $3 \"\t\" $4} ' %s > %s-reformat" %(liftovertempbed, liftovertempbed))
            #    system("mv %s-reformat %s" %(liftovertempbed, liftovertempbed))
            system("bedSort %s %s" %(liftovertempbed, liftovertempbed))
            outbigbed = os.path.join(genomeoutdir, "%s.bb" %othergenome)
            chrsizefile = os.path.join(outdir, othergenome, "chrom.sizes")
            system("bedToBigBed %s %s %s" %(liftovertempbed, chrsizefile, outbigbed))
            #Cleanup:
            system("rm %s" % liftovertempbed)
        system("rm %s" %tempbed) #cleanup

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

def getCGpercent(genomedir, genome):
    twobitfile = os.path.join(genomedir, "%s.2bit" %genome)
    tempfile = os.path.join(genomedir, "%s.gc.wigVarStep.gz" %genome)
    cmd = "hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 %s %s | gzip -c > %s" %(genome, twobitfile, tempfile)
    system(cmd)
    chrsizefile = os.path.join(genomedir, "chrom.sizes")
    gcfile = os.path.join(genomedir, "%s.gc.bw" %genome)
    cmd = "wigToBigWig %s %s %s" %(tempfile, chrsizefile, gcfile)
    system(cmd)
    system("rm -f %s" %tempfile)

def getAlignability(genomedir, genome, halfile):
    outfile = os.path.join(genomedir, "%s.alignability.bw" %genome)
    tempwig = os.path.join(genomedir, "%s.alignability.wig" %genome)
    system("halAlignability %s %s > %s" %(halfile, genome, tempwig))
    chrsizefile = os.path.join(genomedir, "chrom.sizes")
    system("wigToBigWig %s %s %s" %(tempwig, chrsizefile, outfile))
    system("rm -f %s" %tempwig)

def getBasicFiles(halfile, genome2seq2len, outdir):
    for genome, seq2len in genome2seq2len.iteritems():
        genomedir = os.path.join(outdir, genome)
        system("mkdir -p %s" % genomedir)

        makeTwoBitSeqFile(genome, halfile, genomedir) #genomedir/genome.2bit
        getChromSizes(halfile, seq2len, os.path.join(genomedir, "chrom.sizes")) #genomedir/chrom.sizes
        getCGpercent(genomedir, genome) #genomedir/genome.cg.bw
        getAlignability(genomedir, genome, halfile) #genomedir/genome.alignability.bw

def writeTrackDb_cgPercent(f, genome):
    f.write("track cgPercent\n")
    f.write("longLabel GC Percent in 5-base Window\n")
    f.write("shortLabel GC Percent\n")
    f.write("type bigWig 0 100\n")
    f.write("group map\n")
    f.write("visibility dense\n")
    f.write("windowingFunction Mean\n")
    f.write("bigDataUrl %s.gc.bw\n" %genome)
    
    f.write("priority 23.5\n")
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
    
    f.write("priority 23.5\n")
    f.write("autoScale Off\n")
    f.write("maxHeightPixels 128:36:16\n")
    f.write("graphTypeDefault Bar\n")
    f.write("gridDefault OFF\n")
    f.write("color 0,0,0\n")
    f.write("altColor 128,128,128\n")
    f.write("viewLimits 0:%d\n" %genomeCount)
    f.write("\n")
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
    f.write("type bed 3 .")
    f.write("noInherit on\n")
    
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

def writeTrackDb_bigbeds(f, bigbeddir, genomes, currgenome):
    for genome in os.listdir(bigbeddir):
        f.write("track gene%s\n" % genome)
        f.write("longLabel %s Liftovered Genes\n" % genome)
        f.write("shortLabel %sgenes\n" % genome)
        f.write("bigDataUrl ../%s\n" % os.path.join( os.path.basename(bigbeddir), genome, "%s.bb" % currgenome ) )
        f.write("type bigBed 12\n")
        f.write("group annotation\n")
        if genome == currgenome: 
            f.write("visibility dense\n")
        else:
            f.write("visibility hide\n")
        f.write("\n")

def writeTrackDb_snakes(f, halfile, genomes, currgenome):
    for genome in genomes:
        if re.search(genome, currgenome): #current genome
            continue
        #SNAKE TRACKS
        f.write("track snake%s\n" %genome)
        f.write("longLabel %s\n" %genome)
        f.write("shortLabel %s\n" %genome)
        f.write("otherSpecies %s\n" %genome)
        f.write("visibility full\n")
        f.write("bigDataUrl %s\n" % halfile)
        f.write("type halSnake\n")
        f.write("group snake\n")
        f.write("\n")

def writeTrackDbFile(genomes, halfile, outdir, options):

    currgenome = outdir.rstrip('/').split("/")[-1]
    filename = os.path.join(outdir, "trackDb.txt")
    f = open(filename, 'w')
    
    writeTrackDb_cgPercent(f, currgenome)
    writeTrackDb_alignability(f, currgenome, len(genomes))

    if options.bigbeddir:
        writeTrackDb_bigbeds(f, options.bigbeddir, genomes, currgenome)

    if options.rmskdir:
        writeTrackDb_rmsk(f, os.path.join(options.rmskdir, currgenome), outdir)

    writeTrackDb_snakes(f, halfile, genomes, currgenome)
    f.close()

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

def writeGenomesFile(genome2seq2len, halfile, options, outdir):
    '''Write genome for all samples in hal file
    '''

    localHalfile = os.path.join(outdir, os.path.basename(halfile))
    if os.path.abspath(localHalfile) != os.path.abspath(halfile):
        if os.path.exists(localHalfile):
            system("rm %s" %localHalfile)
        if options.cpHal:
            system("cp %s %s" %(os.path.abspath(halfile), localHalfile))
        else:
            system("ln -s %s %s" %(os.path.abspath(halfile), localHalfile))

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
        lodtxtfile, loddir = getLodFiles(localHalfile, options, outdir)
    
    filename = os.path.join(outdir, "genomes.txt")
    f = open(filename, 'w')
    genomes = genome2seq2len.keys()
    for genome in genomes:
        genomedir = os.path.join(outdir, genome)
        f.write("genome %s\n" %genome)
        f.write("twoBitPath %s/%s.2bit\n" % (genome, genome))

        #create trackDb for the current genome:
        if lodtxtfile == '':
            writeTrackDbFile(genomes, "../%s" % os.path.basename(halfile), genomedir, options)
        else:
            writeTrackDbFile(genomes, "../%s" % os.path.basename(lodtxtfile), genomedir, options)
        f.write("trackDb %s/trackDb.txt\n" %genome)
        
        #other info
        f.write("groups groups.txt\n")

        writeDescriptionFile(genome, genomedir)
        f.write("htmlPath %s/description.html\n" %genome)
        f.write("organism %s\n" %genome)
        f.write("orderKey 4800\n")
        f.write("scientificName %s\n" %genome)
        
        seq2len = genome2seq2len[genome]
        (seq, l) = getLongestSeq(seq2len)
        f.write("defaultPos %s:1-%d\n" %(seq, min(l, 1000)))
        f.write("\n")
    f.close()

def writeGroupFile(outdir):
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
    f.write("priority 10\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")
    
    f.write("name annotation\n")
    f.write("label Genes and Other Annotations\n")
    f.write("priority 10\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")
    
    f.write("name exp\n")
    f.write("label Experimental\n")
    f.write("priority 10\n")
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
    parser.add_option('--bedDir', dest='beddir', help='Directory containing bed files of the input genomes. Format: bedDir/ then genome1/ then chr1.bed, chr2.bed... Default=%default' )
    parser.add_option('--rmskDir', dest='rmskdir', help="Directory containing repeatMasker's output files for each genome. Format: rmskDir/ then genome1/ then genome.rmsk.SINE.bb, genome.rmsk.LINE.bb, ... Default=%default")

def checkOptions(parser, args, options):
    if len(args) < 2:
        parser.error("Required two input arguments, %d was provided\n" %len(args))
    if not os.path.exists(args[0]):
        parser.error("Input hal file %s does not exist.\n" %args[0])
    if not os.path.exists(args[1]):
        system("mkdir -p %s" %args[1])
    elif not os.path.isdir(args[1]):
        parser.error("Output directory specified (%s) is not a directory\n" %args[1])

def main():
    usage = '%prog <halFile> <outputDirectory> [options]'
    parser = OptionParser(usage = usage)
    addOptions(parser)
    options, args = parser.parse_args()
    checkOptions(parser, args, options)

    halfile = args[0]
    outdir = args[1]

    writeHubFile(outdir, options)
    writeGroupFile(outdir)
    genomes = getGenomesFromHal(halfile)
    genome2seq2len = getGenomeSequences(halfile, genomes)
  
    #Get basic files (2bit, chrom.sizes, alignability, CG%) for each genome:
    getBasicFiles(halfile, genome2seq2len, outdir)

    #Lift-over gene annotations of each sample to all other samples
    options.bigbeddir = None
    if options.beddir:
        bigbeddir = os.path.join(outdir, "liftoverbeds")
        system("mkdir -p %s" % bigbeddir)
        liftoverBedFiles(options.beddir, halfile, genome2seq2len, bigbeddir, outdir)
        options.bigbeddir = bigbeddir

    writeGenomesFile(genome2seq2len, halfile, options, outdir)

if __name__ == '__main__':
    main()

