#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Mon Oct 28 09:35:55 PDT 2013
Documentation for the conservation track
"""
import sys, os, re, time

def writeConservationDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("This track shows measurements of evolutionary conservation using \n")
    f.write("the <em>phyloP</em> program from the <A HREF=\"http://compgen.bscb.cornell.edu/phast/\" \n")
    f.write("target=_BLANK>PHAST package</A>, for all genomes in the comparative assembly hub.\n")
    f.write("The multiple alignments were generated using \n")
    f.write("<a href=\"https://github.com/glennhickey/progressiveCactus\">progressiveCactus</a>.\n")
    f.write("</P>\n")
    f.write("<P>\n")
    f.write("PhyloP separately measures conservation at individual columns, ignoring the effects of their neighbors.\n") 
    f.write("PhyloP is appropriate for evaluating signatures of selection at particular nucleotides \n")
    f.write("or classes of nucleotides (e.g., third codon positions, or first positions of miRNA target sites).\n")
    f.write("PhyloP can measure acceleration (faster evolution than expected under neutral drift) as well as \n")
    f.write("conservation (slower than expected evolution).  In the phyloP plots, sites \n")
    f.write("predicted to be conserved are assigned positive scores (and shown in blue), \n")
    f.write("while sites predicted to be fast-evolving are assigned negative scores (and \n")
    f.write("shown in red).  The absolute values of the scores represent -log p-values \n")
    f.write("under a null hypothesis of neutral evolution. \n")
    f.write("PhyloP treat alignment gaps and unaligned nucleotides as \n")
    f.write("missing data.\n")
    f.write("</P>\n")
    f.write("\n")

def writeConservationDocs_displays(f):
    f.write("<H2>Display Convention and Configuration</H2>\n")
    f.write("<P>\n")
    f.write("In full and pack display modes, conservation scores are displayed as a \n")
    f.write("<EM>wiggle track</EM> (histogram) in which the height reflects the size of the score.\n")
    f.write("The conservation wiggles can be configured in a variety of ways to\n") 
    f.write("highlight different aspects of the displayed information. \n")
    f.write("Click the <A HREF=\"http://genome.ucsc.edu/goldenPath/help/hgWiggleTrackHelp.html\"\n")
    f.write("TARGET=_blank>Graph configuration help</A> link for an explanation \n") 
    f.write("of the configuration options.</P>\n")
    f.write("\n")

def writeConservationDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("The conservation tracks of this comparative assembly hub were created using the <a href=\"https://github.com/glennhickey/hal/tree/development/phyloP\">phyloP package</a>, which is part of HAL tools.\n")
    f.write("The HAL's <em>phyloP</em> is a python wrapper for running the PHAST package phyloP program and building conservation tracks for all genomes in the HAL multiple alignment.\n")
    f.write("The process starts with creating a neutral model using a reference genome and neutral regions provided.\n")
    f.write("By default, it expects the neutral regions to be coding genes and uses 4fold degenerate site within those genes (which it extracts automatically).\n")
    f.write("The neutral regions could also be ancestral repeats (or anything else).\n")
    f.write("After the model is created, the program proceeds to compute the conservation scores for positions along the root genome.\n")
    f.write("These scores are subsequently lifted over to the children genomes using the multiple alignment.\n")
    f.write("Lastly, HAL phyloP computes the conservation scores for regions in the children genomes that do not align to the root genome.\n")
    f.write("(Of note, the program uses the phylogenetic tree extracted from the HAL file if the tree is not specified.)\n")
    f.write("</P>\n")
    f.write("\n")

def writeConservationDocs_credits(f):
    f.write("<H2>Credits</H2>\n")
    f.write("<P>\n")
    f.write("We thank Melissa Jane Hubisz and Adam Siepel for the PHAST package and their help with HAL phyloP.<br>\n")
    f.write("The HAL phyloP package: Glenn Hickey, Joel Armstrong, Ngan Nguyen, Benedict Paten.<br>\n")
    f.write("</P>\n")
    f.write("\n")

def writeConservationDocs_references(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")
    f.write("Siepel, A., Pollard, K. and Haussler, D.. \n")
    f.write("New methods for detecting lineage-specific selection. \n")
    f.write("<em>Research in Computational Molecular Biology</em>. 2006:190-205.\n")
    f.write("</P>\n")
    f.write("\n")

    f.write("<P>\n")
    f.write("Hickey <i>et al.</i>.\n")
    f.write("<a href=\"http://bioinformatics.oxfordjournals.org/content/29/10/1341.short\">HAL: a hierarchical format for storing and analyzing multiple genome alignments.</a>.\n")
    f.write("<em>Bioinformatics</em>. 2013 May;29(10):1341-1342.\n")
    f.write("</P>\n")
    f.write("\n")

def writeConservationDocs(file):
    f = open(file, 'w')
    writeConservationDocs_description(f)
    writeConservationDocs_displays(f)
    writeConservationDocs_methods(f)
    writeConservationDocs_credits(f)
    writeConservationDocs_references(f)
    f.close()

def makeConservationDocs(outdir):
    outfile = os.path.join(outdir, "conservation.html")
    writeConservationDocs(outfile)

def main():
    makeConservationDocs(sys.argv[1])

if __name__ == '__main__':
    main()

