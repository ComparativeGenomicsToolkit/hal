#!/usr/bin/env python
"""
Generates a single assembly hub that compares two hal files that have genomes in common.

Produces tracks with coverage comparisons between all the genomes, as well as snakes for both hals.

Requires jobTree, sonLib, and bedtools.
"""
import os
from argparse import ArgumentParser

from sonLib.bioio import system, popenCatch, getTempFile
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

def getGenomesInHal(halFile):
    """Get a set of all genomes in the hal file."""
    spaceSeparatedGenomes = popenCatch("halStats --genomes %s" % halFile).strip()
    return set(spaceSeparatedGenomes.split(" "))

def getChromSizes(halFile, genome):
    """Get a dict of sequence name -> sequence size for a particular genome."""
    output = popenCatch("halStats --chromSizes %s %s" % (genome, halFile)).strip()
    splitOutput = [i.split("\t") for i in output.split("\n")]
    return dict((i[0], int(i[1])) for i in splitOutput)

def getGenomeBed(halFile, genome, output):
    """Put a BED file covering all sequences in the genome at the output path."""
    system("halStats --bedSequences %s %s > %s" % (genome, halFile, output))

def createTrackDb(target, genome, genomes, hal1, hal2, label1, label2, hubDir):
    """Create the trackDb.txt for a specific genome."""
    if not os.path.isdir(os.path.join(hubDir, genome)):
        os.makedirs(os.path.join(hubDir, genome))

    # Snake tracks.
    with open(os.path.join(hubDir, genome, 'trackDb.txt'), 'w') as trackDb:
        trackDb.write('''
track alignments
shortLabel Alignments
longLabel Alignments
view Alignments
visibility full
compositeTrack on
type bigBed 3

''')
        for i, targetGenome in enumerate(genomes):
            for halLabel, halPath in zip([label1, label2], [hal1, hal2]):
                trackDb.write('''
	track snake{genome}_{halName}
	longLabel {genome}_{halName}
	shortLabel {genome}_{halName}
	otherSpecies {genome}
	visibility full
	parent alignments
        priority {index}
	bigDataUrl {halPath}
	type halSnake

'''.format(genome=targetGenome, halName=halLabel, index=i, halPath="../" + os.path.basename(halPath)))

def liftoverEntireGenome(target, halFile, fromGenome, toGenome, toBed):
    fromBed = os.path.join(target.getLocalTempDir(), 'from.bed')
    getGenomeBed(halFile, fromGenome, fromBed)
    system("halLiftover %s %s %s %s %s" % (halFile, fromGenome, fromBed, toGenome, toBed))

def subtractAllBeds(target, coverageBeds):
    pass

def subtractBed(target, referenceBed, bedToSubtract, outputPath):
    """Subtract intervals covered by bedToSubtract from referenceBed
    and put the result in outputPath."""
    system("bedtools subtract -a %s -b %s > %s" % (referenceBed, bedToSubtract, outputPath))

def writeHubFile(hubTxtPath, hubName):
    """Write the hub.txt file."""
    with open(hubTxtPath, 'w') as hubFile:
        hubFile.write('''
hub {hubName}
shortLabel {hubName}
longLabel {hubName}
genomesFile genomes.txt
email NoEmail
'''.format(hubName=hubName))

def writeGenomesFile(genomesTxtPath, halFile, genomes):
    """Write the genomes.txt file."""
    with open(genomesTxtPath, 'w') as genomesFile:
        for genome in genomes:
            # Find a valid default chromosome and position. We pick the
            # middle 10000 bases of the maximum-length sequence.
            chromSizes = getChromSizes(halFile, genome)
            maxChrom = max(chromSizes.iterkeys(), key=lambda x: chromSizes[x])
            defaultPosStart = chromSizes[maxChrom] / 2 - 5000
            defaultPosEnd = chromSizes[maxChrom] /  2 + 5000
            defaultPosStart = min(defaultPosStart, 0)
            defaultPosEnd = max(defaultPosEnd, chromSizes[maxChrom])

            genomesFile.write('''
genome {genome}
twoBitPath {genome}/{genome}.2bit
trackDb {genome}/trackDb.txt
organism {genome}
scientificName {genome}
description {genome}
defaultPos {maxChrom}:{defaultPosStart}-{defaultPosEnd}
'''.format(genome=genome, maxChrom=maxChrom, defaultPosStart=defaultPosStart, defaultPosEnd=defaultPosEnd))

def writeSequenceData(genomes, hal, hubDir):
    """Write the .2bit and chrom.sizes for each genome."""
    for genome in genomes:
        if not os.path.isdir(os.path.join(hubDir, genome)):
            os.makedirs(os.path.join(hubDir, genome))

        fasta = getTempFile()
        system("hal2fasta %s %s > %s" % (hal, genome, fasta))
        system("faToTwoBit %s %s" % (fasta, os.path.join(hubDir, genome, genome + '.2bit')))
        system("twoBitInfo %s %s" % (os.path.join(hubDir, genome, genome + '.2bit'), os.path.join(hubDir, genome, 'chrom.sizes')))

def linkHals(hubDir, hal1, hal2):
    """Symlink the hals to the hub directory."""
    system("ln -s %s %s" % (hal1, hubDir))
    system("ln -s %s %s" % (hal2, hubDir))

def createHub(target, genomes, opts):
    """Main method that organizes the creation of the meta-compartive hub."""
    # Create the necessary hub files
    if not os.path.isdir(opts.hubDir):
        os.makedirs(opts.hubDir)
    writeHubFile(os.path.join(opts.hubDir, 'hub.txt'),
                 hubName="%s_vs_%s" % (opts.label1, opts.label2))
    writeGenomesFile(os.path.join(opts.hubDir, 'genomes.txt'), opts.hal1, genomes)
    writeSequenceData(genomes, opts.hal1, opts.hubDir)
    linkHals(opts.hubDir, opts.hal1, opts.hal2)

    # Liftover all genomes
    for genome1 in genomes:
        for genome2 in genomes:
            pass
            # target.addChildTargetFn(liftoverEntireGenome, (opts.hal1, genome1, genome2))
            # target.addChildTargetFn(liftoverEntireGenome, (opts.hal2, genome1, genome2))
    # Create trackDbs
    for genome in genomes:
        target.addChildTargetFn(createTrackDb, (genome, genomes, opts.hal1, opts.hal2, opts.label1, opts.label2, opts.hubDir))
    # Create the bed files that display differential coverage

def parse_args():
    """Parses arguments from sys.argv."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('hal1', help='First hal file')
    parser.add_argument('hal2', help='Second hal file')
    parser.add_argument('hubDir', help='Directory to place the finished hub in')
    parser.add_argument('--label1', help='Label for first hal file (default: hal file name)')
    parser.add_argument('--label2', help='Label for second hal file (default: hal file name)')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()

def main():
    opts = parse_args()
    # Create labels for the HALs if none were provided
    if opts.label1 is None:
        opts.label1 = os.path.basename(opts.hal1)
    if opts.label2 is None:
        opts.label2 = os.path.basename(opts.hal2)

    # Ensure that the hals have some genomes in common
    genomes1 = getGenomesInHal(opts.hal1)
    genomes2 = getGenomesInHal(opts.hal2)
    genomes = genomes1 & genomes2
    if len(genomes) == 0:
        raise ValueError("No genomes in common between the HALs.")

    Stack(Target.makeTargetFn(createHub, (genomes, opts))).startJobTree(opts)

if __name__ == '__main__':
    main()
