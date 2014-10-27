#!/usr/bin/env python
"""Pipeline to get a very simple dotplot relating two sequences in a
HAL file.

Requires R and ggplot2."""
from sonLib.bioio import popenCatch
from argparse import ArgumentParser

def getBedLineForSequence(halFile, genome, sequence):
    """Get a bed line from the beginning to the end of a given
    sequence."""
    bedLines = popenCatch(
        "halStats --bedSequences %s %s" % (genome, halFile)).split("\n")
    seqLines = filter(lambda x: x[0] == sequence, [line.split() for line in bedLines if line != ""])
    if len(seqLines) > 1:
        raise RuntimeError("More than one sequence named %s in genome %s, "
                           "aborting!" % (sequence, genome))
    elif len(seqLines) == 0:
        raise RuntimeError("No sequence named %s found in genome %s" % (sequence,
                                                                        genome))
    return "\t".join(seqLines[0])

def liftoverLine(halFile, refGenome, refBedLine, targetGenome, targetSeq=None):
    """Get a list of PSL lines representing the alignment on the given bed
    line between the refGenome and targetGenome. Optionally, filter
    for only lines involving a certain sequence in targetGenome.
    """
    pslLines = popenCatch("halLiftover --outPSL %s %s stdin %s stdout" % \
                          (halFile,
                           refGenome,
                           targetGenome),
                          stdinString=refBedLine).split("\n")
    pslLines = filter(lambda x: x != "", pslLines)
    if targetSeq is not None:
        pslLines = filter(lambda x: x.split()[13] == targetSeq, pslLines)
    return pslLines

def pslsToDotplotTsv(pslLines, genomeX, seqX, genomeY, seqY):
    """Convert a list of PSL lines to a dotplot TSV comparable to LASTZ's
    dotplot output, namely:
        <target_name>            <query_name>
    <segment1_target_start>  <segment1_query_start>
    <segment1_target_end>    <segment1_query_end>
    NA                       NA
    <segment2_target_start>  <segment2_query_start>
    <segment2_target_end>    <segment2_query_end>
    NA                       NA

    with inclusive 1-based indices, and query end < query start when
    alignments are on different strands.
    """
    tsvLines = ["%s.%s\t%s.%s" % (genomeX, seqX, genomeY, seqY)]
    for line in pslLines:
        fields = line.split()
        numBlocks = int(fields[17])
        (strandX, strandY) = fields[8]
        assert strandX == '+'
        blockSizes = [int(x) for x in fields[18].split(",") if x != ""]
        startsX = [int(x) for x in fields[19].split(",") if x != ""]
        blocksX = [(start, start + size) for start, size in zip(startsX,
                                                                blockSizes)]
        startsY = [int(x) for x in fields[20].split(",") if x != ""]
        blocksY = None
        if strandY == '+':
            blocksY = [(start, start + size) for start, size in zip(startsY,
                                                                    blockSizes)]
        else:
            # have to un-reverse-complement the blocks.
            sizeY = int(fields[14])
            blocksY = [(sizeY - start, sizeY - start - size) for start, size in \
                       zip(startsY, blockSizes)]
        assert len(blocksX) == len(blocksY)
        for blockX, blockY in zip(blocksX, blocksY):
            tsvLines.append("%d\t%d" % (blockX[0], blockY[0]))
            tsvLines.append("%d\t%d" % (blockX[1], blockY[1]))
            tsvLines.append("NA\tNA")
    return tsvLines

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('halFile', help="hal file containing the alignment")
    parser.add_argument('genomeX', help="Genome containing seqX")
    parser.add_argument('seqX', help="Sequence in genomeX to be plotted on X "
                        "axis")
    parser.add_argument('genomeY', help="Genome containing seqY")
    parser.add_argument('seqY', help="Sequence in genomeY to be plotted on Y "
                        "axis")
    opts = parser.parse_args()
    bedLineX = getBedLineForSequence(opts.halFile, opts.genomeX, opts.seqX)
    psls = liftoverLine(opts.halFile, opts.genomeX, bedLineX, opts.genomeY,
                        opts.seqY)
    rInput = pslsToDotplotTsv(psls, opts.genomeX, opts.seqX,
                              opts.genomeY, opts.seqY)
    for line in rInput:
        print line

if __name__ == '__main__':
    main()
