#!/usr/bin/env python
"""
Generates a single assembly hub that compares two hal files that have the same species.

Produces tracks with coverage comparisons between all the genomes, as well as snakes for both hals.

Requires jobTree, sonLib, and bedtools.
"""
from argparse import ArgumentParser

from sonLib.bioio import system, popenCatch, getTempFile
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

def getGenomesInHal(halFile):
    """Get a set of all genomes in the hal file."""
    spaceSeparatedGenomes = popenCatch("halStats --genomes %s" % halFile)
    return set(spaceSeparatedGenomes.split(" "))

def getGenomeBed(halFile, genome, output):
    """Put a BED file covering all sequences in the genome at the output path."""
    system("halStats --bedSequences %s %s > %s" % (genome, halFile, output))

def createTrackDb(target, genome, hubDir):
    """Create the trackDb.txt for a specific genome."""
    pass

def liftoverEntireGenome(target, halFile, fromGenome, toGenome):
    pass

def subtractAllBeds(target, coverageBeds):
    pass

def subtractBed(target, referenceBed, bedToSubtract, outputPath):
    """Subtract intervals covered by bedToSubtract from referenceBed
    and put the result in outputPath."""
    system("bedtools subtract -a %s -b %s > %s" % (referenceBed, bedToSubtract, outputPath))

def createHub(target, genomes, opts):
    """Main method that organizes the creation of the meta-compartive hub."""
    # Liftover all genomes
    for genome1 in genomes:
        for genome2 in genomes:
            target.addChildTargetFn(liftoverEntireGenome, (opts.hal1, genome1, genome2))
            target.addChildTargetFn(liftoverEntireGenome, (opts.hal2, genome1, genome2))
    # Create trackDbs
    for genome in genomes:
        target.addChildTargetFn(createTrackDb, (genome, opts.hubDir))
    # Create the bed files that display differential coverage
    target.setFollowOnTargetFn(subtractAllBeds, (coverageBed))

def parse_args():
    """Parses arguments from sys.argv."""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('hal1', help='First hal file')
    parser.add_argument('hal2', help='Second hal file')
    parser.add_argument('hubDir', help='Directory to place the finished hub in')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()

def main():
    opts = parse_args()

    # Ensure that the hals cover the same genomes
    genomes1 = getGenomesFromHal(opts.hal1)
    genomes2 = getGenomesFromHal(opts.hal2)
    if genomes1 != genomes2:
        raise ValueError('Hal files do not contain the same genomes.')
    genomes = genomes1

    Stack(Target.makeTargetFn(createHub, (genomes, opts))).startJobTree(opts)

if __name__ == '__main__':
    main()
