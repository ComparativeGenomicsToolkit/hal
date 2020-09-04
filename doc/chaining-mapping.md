# Pairwise chaining and alignment mapping with HAL alignments

This document describes how to used produce UCSC browser-style pairwise chains 
from HAL alignments.  It also describes how to use these chains project alignments
and annotations between genomes (we refer to this as TransMap).

\subsection{Human/chicken transcript alignment protocols}\label{sec:tanscriptMappingProtocol}
Protein-coding transcript annotations were obtained from the UCSC Genome browser
\cite{ucscBrowser} tables.  Human annotations are GENCODE V34 on hg38
(GRCh38/GCA\_000001405.27) and chicken annotations are Ensembl 85 on galGal4
(GCA\_000002315.2).  Predicted RNA sequences for each protein-coding transcript are extracted
from the genome.  Only gene annotations on the primary assemblies were used,
those on alternate loci, patches, and assembled sequences were dropped.
This results in 84,001 transcripts in 19,695 genes for human and
15,328 transcript in 14,499 genes for chicken.
These are then mapped from the source genome to the destination genome.

The steps for each method are outlined below, although the actual execution was
done by partitioning the data and using a cluster. Command-line tools from the
UCSC Genome Browser group are used for processing the resulting alignments.

\subsubsection{BLATX transcript alignment protocol}
The \textit{BLATX} alignments were created using protein-translated
mode to align the mRNAs to the target genome with \textit{BLAT}~\cite{BLAT}
version 36x5.  They were then filtered following the same protocol the UCSC
Genome Browser uses for creating the other species RefSeq alignments:

\begin{verbatim}
    blat -noHead -q=rnax -t=dnax -mask=lower <dest-genome.2bit> \
        <src-rna.fa> <dest-rna-raw.psl>
\end{verbatim}

We then filter to get near-best in genome.  Alignment to chicken uses near best filter
of \code{-localNearBest=0.010} while to human it uses \code{-globalNearBest=0.010}:

\begin{verbatim}
    faPolyASizes <src-rna.fa> <src-rna.polya>
    pslCDnaFilter <nearBestOption> -minId=0.35 -minCover=0.15 -minQSize=20 \
        -ignoreIntrons -repsAsMatch -ignoreNs -bestOverlap \
        -polyASizes=<src-rna.polya> <dest-rna-raw.psl> <dest-rna-mapped.psl>
\end{verbatim}

The \code{transMapPslToGenePred} command is then used to project the original CDS
onto the alignment.

\subsubsection{TBLASTX transcript alignment protocol}
The \textit{TBLASTX} alignments were created using protein-translated
\textit{tblastx} program to align the mRNAs to the target genome with
\textit{BLAST+}~\cite{BLASTPLUS} version 2.10.0+.

The database is created using the repeat masking from the UCSC Genome Browser
genomes to match what is used with in the \textit{BLATX} methodology above:
\begin{verbatim}
    convert2blastmask -in <dest-genome.fa> -masking_algorithm repeat \
        -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin \
        -out <dest-genome.mask>
    makeblastdb -dbtype nucl -in <dest-genome.fa> -mask_data <dest-genome.mask>
\end{verbatim}

The mRNAs are aligned and the resulting XML converted to PSL format, filtering to an
e-value threshold of 0.01.  These are then chained using a program the UCSC group
developed for chaining \textit{BLAST} alignments:
\begin{verbatim}
    tblastx -db <dest-genome.fa> -db_soft_mask 40 -outfmt 5 -query <src-rna.fa> \
        -out <dest-rna-raw.xml>
    blastXmlToPsl -eVal=0.01 <dest-rna-raw.xml> <dest-rna-raw.psl>
    simpleChain -outPsl -maxGap=75000 <dest-rna-raw.psl> <dest-rna-chained.psl>
\end{verbatim}

The alignments produced are then filtered in the same manner as the \textit{BLATX}
alignments.


\subsubsection{LASTZ transcript alignment protocol}
Both the \textit{LASTZ}~\cite{lastz} and \textit{Cactus} transcript mappings uses the
\textit{TransMap}~\cite{transMap} projection alignment algorithm to project
transcript annotation between genomes.

The \textit{LASTZ} alignment chains and nets
\cite{CHIAROMONTE_2001,lastz,kent2003evolution,Schwartz_2003} were obtained from the
UCSC Genome Browser downloads.  These were then filtered to produce a set of
syntenic mapping chains using these steps:

\begin{verbatim}
    netFilter -syn <genomes.net> <syntenic.net>
    netChainSubset -wholeChains <syntenic.net> <genome.chain> <mapping.chain>
\end{verbatim}

\subsubsection{Cactus transcript alignment protocol}
The \textit{Cactus} alignments are extracted for all primary chromosomes from the HAL
file and chained using the same chaining algorithm as the \textit{LASTZ} chains, with
the \code{--noDupes} option having a similar effect as the syntenic net filtering:
\begin{verbatim}
   halLiftover --outPSL --noDupes 600way.hal <srcOrganism> \
      <srcChroms.bed> <destOrganism>  <src-dest.psl> <genome.psl>
   axtChain -psl -linearGap=loose -scoreScheme=HoxD55.q <genome.psl> <mapping.chain>
\end{verbatim}

The \textit{TransMap} protocol is use for both the \textit{LASTZ} and
\textit{Cactus} mapping chains to produce alignments of the transcripts to the
other genomes.  This used the \textit{pslMap} command to do the mapping and
\textit{pslRecalcMatch} to update the statistic in the alignments:

\begin{verbatim}
    pslMap -chainMapFile <src-rna.psl> <mapping.chains> <dest-rna-over.psl>
    pslRecalcMatch <dest-rna-over.psl> <dest-genome.2bit> <src-rna.fa> <dest-rna-raw.psl>
\end{verbatim}

The alignments produced are then filtered in the same manner as the \textit{LASTZ}
alignments.

