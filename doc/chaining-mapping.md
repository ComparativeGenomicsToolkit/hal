# FIXME: remove --noDupes

# Pairwise chaining and alignment mapping with HAL alignments

This document describes how to produce UCSC browser-style pairwise chains from HAL alignments.  It also describes how to use these chains to project alignments or annotations between genomes (we refer to this as TransMap).

##  Creating chains from HAL

Creating chains uses the `halLiftover` program to produce pairwise alignments for each chromosome in PSL format.  These alignments are then chained using the UCSC browser chaining protocol.

The names of the source (query) and target genomes are required, as specified in the HAL. The list of genome names in a HAL command.

```
halStats --genomes 600way.hal
```

The source and target genome sequences in UCSC two-bit format are required and obtained with:

```
hal2fasta 600way.hal Homo_sapiens | faToTwoBit stdin Homo_sapiens.2bit
hal2fasta 600way.hal Gallus_gallus | faToTwoBit stdin Gallus_gallus.2bit
```

Next, the source genome sequence name full ranges are  obtained in BED format:

```
halStats --bedSequences Homo_sapiens 600way.hal > Homo_sapiens.bed
```

The BED is then used by ``halLiftover`` to create pairwise alignments,
which are forced to the positive strand with `pslPosTarget`:

```
halLiftover --outPSL --noDupes 600way.hal Homo_sapiens \
      Homo_sapiens.bed Gallus_gallus /dev/stdout | \
      pslPosTarget stdin Homo_sapiens-Gallus_gallus.psl
```

The `halLiftover` may be run on a cluster by splitting the source BED into one sequence per job.  The resulting PSL files are then concatenated together.

These alignments are then chained using the UCSC Browser `axtChain` program:

```
axtChain -psl -linearGap=loose Homo_sapiens-Gallus_gallus.psl Gallus_gallus.2bit Homo_sapiens.2bit Homo_sapiens-Gallus_gallus.chain
```

Note the order of the two-bit files are in the order target, followed by source.  Consult the output of `axtChain` with no arguments for a description of the `-linearGap` option. It may also improve results to use the `-scoreScheme` option, although find the score files to is challenging.

# Projection alignments using the TransMap protocol.

coming soon
