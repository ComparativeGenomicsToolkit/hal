#!/usr/bin/env python3
"""
Generates a single assembly hub that compares two hal files that have genomes in common.

Produces tracks with coverage comparisons between all the genomes, as well as snakes for both hals.

Requires jobTree, sonLib, and bedtools.
"""
import os
from argparse import ArgumentParser

from sonLib.bioio import system, popenCatch, getTempFile
from toil.job import Job
from toil.common import Toil
from functools import reduce

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

def createTrackDb(target, genome, genomes, hals, labels, hubDir):
    """Create the trackDb.txt for a specific genome."""
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
subGroup1 view Track_Type Snake=Alignments
subGroup2 orgs Organisms {genomeSubgroupList}
dimensions dimensionX=view dimensionY=orgs

'''.format(genomeSubgroupList=" ".join([genome + "=" + genome for genome in genomes])))
        for i, targetGenome in enumerate(genomes):
            for halLabel, halPath in zip(labels, hals):
                trackDb.write('''
	track snake{genome}_{halName}
	longLabel {genome}_{halName}
	shortLabel {genome}_{halName}
	otherSpecies {genome}
	visibility full
	parent alignments
        priority {index}
        subGroups view=Snake orgs={genome}
	bigDataUrl {halPath}
	type halSnake

'''.format(genome=targetGenome, halName=halLabel, index=i, halPath=halPath))

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
            maxChrom = max(iter(chromSizes.keys()), key=lambda x: chromSizes[x])
            defaultPosStart = chromSizes[maxChrom] / 2 - 5000
            defaultPosEnd = chromSizes[maxChrom] /  2 + 5000
            defaultPosStart = max(defaultPosStart, 0)
            defaultPosEnd = min(defaultPosEnd, chromSizes[maxChrom])

            genomesFile.write('''
genome {genome}
twoBitPath {genome}/{genome}.2bit
trackDb {genome}/trackDb.txt
organism {genome}
scientificName {genome}
description {genome}
defaultPos {maxChrom}:{defaultPosStart}-{defaultPosEnd}
'''.format(genome=genome, maxChrom=maxChrom, defaultPosStart=defaultPosStart, defaultPosEnd=defaultPosEnd))

def writeSequenceData(target, genome, hal, hubDir):
    """Write the .2bit and chrom.sizes for a genome."""
    if not os.path.isdir(os.path.join(hubDir, genome)):
        os.makedirs(os.path.join(hubDir, genome))
    fasta = getTempFile()
    system("hal2fasta %s %s > %s" % (hal, genome, fasta))
    system("faToTwoBit %s %s" % (fasta, os.path.join(hubDir, genome, genome + '.2bit')))
    system("twoBitInfo %s %s" % (os.path.join(hubDir, genome, genome + '.2bit'), os.path.join(hubDir, genome, 'chrom.sizes')))
    os.remove(fasta)

def linkHals(hubDir, hals):
    """Symlink the hals to the hub directory."""
    relativePaths = []
    for i, hal in enumerate(hals):
        uniqueName = 'input_%d' % i + ".hal"
        system("ln -sf %s %s" % (hal, os.path.join(hubDir, uniqueName)))
        relativePaths.append('../' + uniqueName)
    return relativePaths

def createHub(target, genomes, opts):
    """Main method that organizes the creation of the meta-compartive hub."""
    # Create the necessary hub files
    if not os.path.isdir(opts.hubDir):
        os.makedirs(opts.hubDir)
    writeHubFile(os.path.join(opts.hubDir, 'hub.txt'),
                 hubName="_vs_".join(opts.labels))
    writeGenomesFile(os.path.join(opts.hubDir, 'genomes.txt'), opts.hals[0], genomes)
    for genome in genomes:
        target.addChildFn(writeSequenceData, (genome, opts.hals[0], opts.hubDir))
    relativeHalPaths = linkHals(opts.hubDir, opts.hals)

    # Liftover all genomes
    for genome1 in genomes:
        for genome2 in genomes:
            pass
            # target.addChildFn(liftoverEntireGenome, (opts.hal1, genome1, genome2))
            # target.addChildFn(liftoverEntireGenome, (opts.hal2, genome1, genome2))
    # Create trackDbs
    for genome in genomes:
        target.addChildFn(createTrackDb, (genome, genomes, relativeHalPaths, opts.labels, opts.hubDir))
    # Create the bed files that display differential coverage

def parse_args():
    """Parses arguments from sys.argv."""
    parser = ArgumentParser(description=__doc__)
    Job.Runner.addToilOptions(parser)
    parser.add_argument('hubDir', help='Directory to place the finished hub in')
    parser.add_argument('hals', type=os.path.abspath, nargs='+', help='Hal files')
    parser.add_argument('--labels', nargs='+', help='Labels for hal files (default: hal file name)')

    return parser.parse_args()

def main():

    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    opts = parse_args()
    opts.hubDir = os.path.abspath(opts.hubDir)
    opts.hals = [os.path.abspath(hal) for hal in opts.hals]
    if opts.batchSystem != 'singleMachine':
        raise RuntimeError("singleMachine is the only supported batchSystem")

    # Create labels for the HALs if none were provided
    if opts.labels is None:
        opts.labels = [os.path.basename(hal) for hal in opts.hals]
    if len(opts.labels) != len(opts.hals):
        raise ValueError("%d labels were provided, but %d hals were provided." % (len(opts.labels), len(opts.hals)))

    # Ensure that the hals have some genomes in common, and take the
    # common genomes to display in the hub.
    genomess = [getGenomesInHal(hal) for hal in opts.hals]
    genomes = reduce(lambda a, i: a.intersection(i), genomess)
    if len(genomes) == 0:
        raise ValueError("No genomes in common between the HALs.")

    with Toil(opts) as toil:
        toil.start(Job.wrapJobFn(createHub, genomes, opts))
    
if __name__ == '__main__':
    main()
