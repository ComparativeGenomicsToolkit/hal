#!/usr/bin/env python3
"""Pipeline to get a very simple dotplot relating two sequences in a
HAL file.

Requires R and ggplot2."""
from sonLib.bioio import popenCatch
from argparse import ArgumentParser

def in_range(element, r):
    """because python's xrange.__contains__ is O(n) for no good reason.
    Doesn't work if element is "out of phase" with r's step.
    """
    return element >= r[0] and element <= r[-1]

def getBedLineForSequence(halFile, genome, sequence, start, length):
    """Get a bed line from the beginning to the end of a given
    sequence. If start and length are None, the full sequence is
    returned, otherwise only the given region is.
    """
    bedLines = popenCatch(
        "halStats --bedSequences %s %s" % (genome, halFile)).split("\n")
    seqLines = [x for x in [line.split() for line in bedLines if line != ""] if x[0] == sequence]
    if len(seqLines) > 1:
        raise RuntimeError("More than one sequence named %s in genome %s, "
                           "aborting!" % (sequence, genome))
    elif len(seqLines) == 0:
        raise RuntimeError("No sequence named %s found in genome %s" % (sequence,
                                                                        genome))
    if start is None and length is None:
        return "\t".join(seqLines[0])
    elif start is not None and length is not None:
        if start + length > int(seqLines[0][2]):
            raise RuntimeError("Selected region runs off end of sequence.")
        seqLines[0][1] = str(start)
        seqLines[0][2] = str(start + length)
        return "\t".join(seqLines[0])
    else:
        raise RuntimeError("Both start and length must be provided.")

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
    pslLines = [x for x in pslLines if x != ""]
    if targetSeq is not None:
        pslLines = [x for x in pslLines if x.split()[13] == targetSeq]
    return pslLines

def pslsToDotplotTsv(pslLines, genomeX, seqX, genomeY, seqY, startY, lengthY):
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

        if startY is not None and lengthY is not None:
            # Filter out blocks that are outside the selected region in Y
            # (the filtering for the selected region in X took place
            # before the liftover).
            regionY = range(startY, startY + lengthY)
            deleteIndices = [index for index, block in enumerate(blocksY) if not in_range(block[0], regionY) or not in_range(block[1], regionY)]
            numDeleted = 0
            for deleteIndex in deleteIndices:
                del blocksY[deleteIndex - numDeleted]
                del blocksX[deleteIndex - numDeleted]
                numDeleted += 1

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
    parser.add_argument('--startX', help="Start position in sequence X",
                        type=int)
    parser.add_argument('--startY', help="Start position in sequence Y",
                        type=int)
    parser.add_argument('--lengthX', help="Length in sequence X",
                        type=int)
    parser.add_argument('--lengthY', help="Length in sequence Y",
                        type=int)
    opts = parser.parse_args()
    bedLineX = getBedLineForSequence(opts.halFile, opts.genomeX, opts.seqX,
                                     opts.startX, opts.lengthX)
    psls = liftoverLine(opts.halFile, opts.genomeX, bedLineX, opts.genomeY,
                        opts.seqY)
    rInput = pslsToDotplotTsv(psls, opts.genomeX, opts.seqX,
                              opts.genomeY, opts.seqY, opts.startY,
                              opts.lengthY)
    for line in rInput:
        print(line)

if __name__ == '__main__':
    main()
