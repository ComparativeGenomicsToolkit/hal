#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Mon Oct 28 09:35:55 PDT 2013
Documentation for the repeatMasker track
"""
import sys, os

def writeRepeatMaskerDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("This track was created by using Arian Smit's <a href=\"http://www.repeatmasker.org/\">RepeatMasker</a> program, which screens DNA sequences for interspersed repeats and low complexity DNA sequences. The program outputs a detailed annotation of the repeats that are present in the query sequence (represented by this track), as well as a modified version of the query sequence in which all the annotated repeats have been masked. RepeatMasker uses the <a href=\"http://www.girinst.org/repbase/update/index.html\">Repbase Update</a> library of repeats from the <a href=\"http://www.girinst.org/\">Genetic Information Research Institute</a> (GIRI). Repbase Update is described in Jurka (2000) in the References section below.\n")
    f.write("</P>\n")
    f.write("\n")

def writeRepeatMaskerDocs_displays(f):
    f.write("<H2>Display Conventions and Configuration</H2>\n")
    f.write("<P>\n")
    f.write("In full display mode, this track displays up to ten different classes of repeats:\n")
    f.write("<ul>\n")
    f.write("<li>Short interspersed nuclear elements (SINE), which include ALUs</li>\n")
    f.write("<li>Long interspersed nuclear elements (LINE)</li>\n")
    f.write("<li>Long terminal repeat elements (LTR), which include retroposons</li>\n")
    f.write("<li>DNA repeat elements (DNA)</li>\n")
    f.write("<li>Simple repeats (micro-satellites)</li>\n")
    f.write("<li>Low complexity repeats</li>\n")
    f.write("<li>Satellite repeats</li>\n")
    f.write("<li>RNA repeats (including RNA, tRNA, rRNA, snRNA, scRNA, srpRNA)</li>\n")
    f.write("<li>Other repeats, which includes class RC (Rolling Circle)</li>\n")
    f.write("<li>Unknown</li>\n")
    f.write("</ul>\n")
    f.write("</P>\n")

    f.write("<p>\n")
    f.write("The level of color shading in the graphical display reflects the amount of \
    base mismatch, base deletion, and base insertion associated with a repeat element. \
    The higher the combined number of these, the lighter the shading.\n")
    f.write("</p>\n")

    f.write("<p>\n")
    #f.write("A &quot;?&quot; at the end of the &quot;Family&quot; or &quot;Class&quot; \
    #(for example, DNA?) signifies that the curator was unsure of the classification. \
    #At some point in the future, either the &quot;?&quot; will be removed or the \
    #classification will be changed.\n")
    f.write("</p>\n")
    f.write("\n")

def writeRepeatMaskerDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("Data are generated using the RepeatMasker. Repeats are soft-masked. \
    Alignments may extend through repeats, but are not permitted to initiate in them.\n")
    f.write("</P>\n")
    f.write("\n")

def writeRepeatMaskerDocs_credits(f):
    f.write("<H2>Credits</H2>\n")
    f.write("<P>\n")
    f.write("Thanks to Arian Smit, Robert Hubley and GIRI for providing the tools and\
    repeat libraries used to generate this track.\n")
    f.write("</P>\n")
    f.write("\n")

def writeRepeatMaskerDocs_references(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")

    f.write("<p>\n")
    f.write("Smit AFA, Hubley R, Green P. <em>RepeatMasker Open-3.0</em>.\n")
    f.write("<a href=\"http://www.repeatmasker.org\" target=\"_blank\">\n")
    f.write("http://www.repeatmasker.org</a>. 1996-2010.\n")
    f.write("</p>\n")

    f.write("<p>\n")
    f.write("Repbase Update is described in:\n")
    f.write("</p>\n")

    f.write("<p>\n")
    f.write("Jurka J.\n")
    f.write("<a href=\"http://www.sciencedirect.com/science/article/pii/S016895250002093X\" target=\"_blank\">\n")
    f.write("Repbase Update: a database and an electronic journal of repetitive elements</a>.\n")
    f.write("<em>Trends Genet</em>. 2000 Sep;16(9):418-420.\n")
    f.write("PMID: <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/10973072\" target=\"_blank\">10973072</a>\n")
    f.write("</p>\n")

    f.write("<p>\n")
    f.write("For a discussion of repeats in mammalian genomes, see:\n")
    f.write("</p>\n")

    f.write("<p>\n")
    f.write("Smit AF.\n")
    f.write("<a href=\"http://www.sciencedirect.com/science/article/pii/S0959437X99000313\" target=\"_blank\">\n")
    f.write("Interspersed repeats and other mementos of transposable elements in mammalian genomes</a>.\n")
    f.write("<em>Curr Opin Genet Dev</em>. 1999 Dec;9(6):657-63.\n")
    f.write("PMID: <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/10607616\" target=\"_blank\">10607616</a>\n")
    f.write("</p>\n")

    f.write("<p>\n")
    f.write("Smit AF.\n")
    f.write("<a href=\"http://www.sciencedirect.com/science/article/pii/S0959437X9680030X\" target=\"_blank\">\n")
    f.write("The origin of interspersed repeats in the human genome</a>.\n")
    f.write("<em>Curr Opin Genet Dev</em>. 1996 Dec;6(6):743-8.\n")
    f.write("PMID: <a href=\"http://www.ncbi.nlm.nih.gov/pubmed/8994846\" target=\"_blank\">8994846</a>\n")
    f.write("</p>\n")

    f.write("\n")

def writeRepeatMaskerDocs(file):
    f = open(file, 'w')
    writeRepeatMaskerDocs_description(f)
    writeRepeatMaskerDocs_displays(f)
    writeRepeatMaskerDocs_methods(f)
    writeRepeatMaskerDocs_credits(f)
    writeRepeatMaskerDocs_references(f)
    f.close()

def makeRepeatMaskerDocs(outdir):
    outfile = os.path.join(outdir, "repeatMasker.html")
    writeRepeatMaskerDocs(outfile)

def main():
    makeRepeatMaskerDocs(sys.argv[1])

if __name__ == '__main__':
    main()

