#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Mon Oct 28 09:35:55 PDT 2013
Documentation for the gcPercent track
"""
import sys, os

def writeGcPercentDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("The GC percent track shows the percentage of G (guanine) and C (cytosine) bases in 5-base windows. High GC content is typically associated with gene-rich areas.\n")
    f.write("</P>\n")
    f.write("<P>\n")
    f.write("This track may be configured in a variety of ways to highlight different aspects of the displayed information. Click the \"Graph configuration help\" link for an explanation of the configuration options.\n")
    f.write("</P>\n")
    f.write("\n")

def writeGcPercentDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("This track was generated following the <a href=\"http://genomewiki.ucsc.edu/index.php/Browser_Track_Construction#GC_Percent\">UCSC GC_Percent Track Construction instructions</a>, using the sequence information extracted from the multiple sequence alignments.\n")
    f.write("</P>\n")
    f.write("\n")

def writeGcPercentDocs_credits(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")
    f.write("The GC Percent graph presentation is by <a href=\"mailto:hiram@soe.ucsc.edu\">Hiram Clawson</a>. The data was automatically generated using the <a href=\"https://github.com/glennhickey/hal\">HAL tools package</a>.\n")
    f.write("</P>\n")
    f.write("\n")

def writeGcPercentDocs(file):
    f = open(file, 'w')
    writeGcPercentDocs_description(f)
    writeGcPercentDocs_methods(f)
    writeGcPercentDocs_credits(f)
    f.close()

def makeGcPercentDocs(outdir):
    outfile = os.path.join(outdir, "gcPercent.html")
    writeGcPercentDocs(outfile)

def main():
    makeGcPercentDocs(sys.argv[1])

if __name__ == '__main__':
    main()

