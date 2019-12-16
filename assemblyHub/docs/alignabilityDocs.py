#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Mon Oct 28 09:35:55 PDT 2013
Documentation for the alignability track
"""
import sys, os

def writeAlignabilityDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("This track shows the number of genomes aligned to each position of the reference. The values range from 0 to the total number of input genomes and imputed ancestral genomes.\n")
    f.write("</P>\n")
    f.write("\n")

def writeAlignabilityDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("Alignability was generated using the <em>halAlignability</em> script of the <a href=\"https://github.com/glennhickey/hal\">HAL tools package</a>.\n") 
    f.write("</P>\n")
    f.write("\n")

def writeAlignabilityDocs_references(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")
    f.write("Hickey <i>et al.</i>.\n")
    f.write("<a href=\"http://bioinformatics.oxfordjournals.org/content/29/10/1341.short\">HAL: a hierarchical format for storing and analyzing multiple genome alignments.</a>.\n")
    f.write("<em>Bioinformatics</em>. 2013 May;29(10):1341-1342.\n")
    f.write("</P>\n")
    f.write("\n")

def writeAlignabilityDocs(file):
    f = open(file, 'w')
    writeAlignabilityDocs_description(f)
    writeAlignabilityDocs_methods(f)
    writeAlignabilityDocs_references(f)
    f.close()

def makeAlignabilityDocs(outdir):
    outfile = os.path.join(outdir, "alignability.html")
    writeAlignabilityDocs(outfile)

def main():
    makeAlignabilityDocs(sys.argv[1])

if __name__ == '__main__':
    main()


