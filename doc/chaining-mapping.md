# Pairwise chaining and alignment mapping with HAL alignments

This document describes how to produce UCSC browser-style pairwise chains from
HAL alignments.  It also describes how to use these chains to project
alignments or annotations between genomes.  Chaining outside of HAL is necessary,
as the pairwise alignments extract from HAL are fragmented and lack long-range
continuity.

The goal is to generate UCSC alignment chains to use in projecting alignments
or annotations on a *source* genome to a *destination* genome.   With a genomic
alignment of genome *A* to genome *B*, one can project and annotation or 
alignment of on genome *A* of *X*<sub>*A*</sub>. to genome *B*, producing *X*<sub>*B*</sub>.

##  Creating chains from HAL

UCSC pairwise alignment formats (PSL and chain) alignments have the two sides
of the alignment are the *query* and *target*.  Which of these is consider the
*source* and *destination* alignment differs depending on the tool that is used.
The UCSC `liftOver` program considered the *target* to be the source, so the these
instructions create alignments where the *target* is the *source* genome.

Creating chains uses the `halLiftover` program to produce pairwise alignments
for each chromosome in PSL format.  These alignments are then chained using
the UCSC browser chaining protocol.

The names of the *source* and *destination* genomes are required, as specified
in the HAL. The list of genome names in a HAL command.

```
halStats --genomes 600way.hal
```

The *source* and *destination* genome sequences in UCSC two-bit format are
required and obtained with:

```
hal2fasta 600way.hal Homo_sapiens | faToTwoBit stdin Homo_sapiens.2bit
hal2fasta 600way.hal Gallus_gallus | faToTwoBit stdin Gallus_gallus.2bit
```

Next, the *source* genome sequences are obtained in BED format:

```
halStats --bedSequences Homo_sapiens 600way.hal > Homo_sapiens.bed
```

The BED is then used by ``halLiftover`` to create pairwise alignments,
which are forced to the positive strand with `pslPosTarget`:

```
halLiftover --outPSL 600way.hal Homo_sapiens \
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
