#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Mon Oct 28 09:35:55 PDT 2013
Documentation for the hubCentral composite track
"""
import sys, os, re, time

def writeSnakeDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("An alignment track, or <em>snake</em> track, shows the relationship between the chosen browser genome, termed the reference (genome), and another genome, termed the query (genome).\n")
    f.write("The <em>snake</em> display is capable of showing all possible types of structural rearrangement.\n")
    f.write("</P>\n")
    f.write("\n")

def writeSnakeDocs_displays(f):
    f.write("<H2>Display Convention and Configuration</H2>\n")
    f.write("<P>\n")
    f.write("In <em>full</em> display mode, a <em>snake</em> track can be decomposed into two primitive drawing elements, segments, which are the colored rectangles, and adjacencies, which are the lines connecting the segments. Segments represent subsequences of the query genome aligned to the given portion of the reference genome. Adjacencies represent the covalent bonds between the aligned subsequences of the query genome. Segments can be configured to be colored by chromosome, strand or left a single color under the <em>Select track Type</em>, <em>Alignments</em>, then <em>Block coloring method</em>.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("Red tick-marks within segments represent substitutions with respect to the reference, shown in windows of the reference of (by default) up to 50 kilo-bases. This default can be adjusted under <em>Select track Type</em>, <em>Alignments</em>, then <em>Maximum window size in which to show mismatches</em>. Zoomed in to the base-level these substitutions are labeled with the non-reference base.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("An insertion in the reference relative to the query creates a gap between abutting segment sides that is connected by an adjacency. An insertion in the query relative to the reference is represented by an orange tick mark that splits a segment at the location the extra bases would be inserted. Simultaneous independent insertions in both query and reference look like an insertion in the reference relative to the query, except that the corresponding adjacency connecting the two segments is colored orange. More complex structural rearrangements create adjacencies that connect the sides of non-abutting segments in a natural fashion.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("Duplications within the query genome create extra segments that overlap along the reference genome axis. Duplications within the reference imply self-alignments, intervals of the reference genome that align to other intervals of the reference genome. To show these self-alignments within the reference genome we draw colored coded sets of lines along the reference genome axis that indicate these self homologies, and align any query segments that align to these regions arbitrarily to just one copy of the reference self alignment.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("The <em>pack</em> display option can be used to display a larger number of Snake tracks in limited vertical browser. This mode eliminates the adjacencies from the display and forces the segments onto as few rows as possible, given the constraint of still showing duplications in the query sequence.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("The <em>dense</em> display further eliminates these duplications so that each Snake track is compactly represented along just one row.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("To ensure that the <em>snake</em> alignments track loads quickly at any resolution, from windows showing individual bases up to entire scaffolds or chromosomes, the LOD (Levels-Of-Detail) algorithm (part of the HAL tools package) is used, which creates scaleable levels of detail for the alignments. The additional use of the hdf5 caching scheme further aides scaling.\n")
    f.write("</P>\n")
    
    f.write("<P>\n")
    f.write("Various mouse overs are implemented and clicking on segments navigates to the corresponding region in the query genome, making it simple to instantly switch the alignment view between reference points.\n")
    f.write("</P>\n")
    
    f.write("\n")

def writeSnakeDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("A <em>snake</em> is a way of viewing a set of pairwise gap-less alignments that may overlap on both the reference and query genomes. Alignments are always represented as being on the positive strand of the reference species, but can be on either strand on the query sequence.\n")
    f.write("</P>\n")
    f.write("\n")
    f.write("<P>\n")
    f.write("A <em>snake</em> plot puts all the query segments within a reference chromosome range on a set of one or more levels. All the segments on a level are on the same strand, do not overlap in reference coordinate space, and are in the same order and orientation in both sequences. This is the same requirement as the alignments in a chain on the UCSC browser. Before the algorithm is started, all the segments are sorted by their starting coordinate on the query, and the current level is set to one. Then in a recursive fashion, the algorithm places the first segment on the current list on the current level, and then adds all the rest of the segments on the list that will fit onto the current level with the requirements that all the segments on a level are on the same strand, and that the proposed segment be non-overlapping and have a reference start address that is greater than the query end address of the previously added segment on that level. All segments that will not fit on the current level are then added to subsequent levels following the same rules. Once all the segments have been assigned a level, lines are drawn between the segments to show the adjacencies in the list when sorted by query start address.\n")
    f.write("</P>\n")
    f.write("\n")

def writeSnakeDocs_credits(f):
    f.write("<H2>Credits</H2>\n")
    f.write("<P>\n")
    f.write("The <em>snake</em> alignment display was implemented by <a href=\"mailto:braney@soe.ucsc.edu\">Brian Raney</a>.<br>\n")
    f.write("HAL supports and track generations: Glenn Hickey, Ngan Nguyen, Joel Armstrong, Benedict Paten.\n")
    f.write("</P>\n")
    f.write("\n")

def writeSnakeDocs_references(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")
    f.write("Paten <i>et al.</i>.\n")
    f.write("<a href=\"http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166836/\">Cactus: Algorithms for genome multiple sequence alignment.</a>\n")
    f.write("<em>Genome research</em>. 2011;21:1512-1528.\n")
    f.write("</P>\n")
    f.write("\n")

def writeSnakeDocs(f):
    f.write("<H1>Alignments</H1>\n")
    writeSnakeDocs_description(f)
    writeSnakeDocs_displays(f)
    writeSnakeDocs_methods(f)
    writeSnakeDocs_credits(f)
    writeSnakeDocs_references(f)

def writeLiftoverDocs_description(f):
    f.write("<H2>Description</H2>\n")
    f.write("<P>\n")
    f.write("Lifted-over annotation tracks show the annotations of any genome translated onto the reference genome, via a process of lift-over.\n")
    f.write("All the alignments and lifted over annotations shown are mutually consistent with one another, because the annotation lift over and alignment display is symmetrically driven by one reference free alignment process, rather than a mixture of different pairwise and reference based multiple alignments.\n")
    f.write("</P>\n")
    f.write("\n")

def writeLiftoverDocs_methods(f):
    f.write("<H2>Methods</H2>\n")
    f.write("<P>\n")
    f.write("The lifted-over tracks were generated using the <em>halLiftover</em> and/or the <em>halWiggleLiftover</em> scripts of the <a href=\"https://github.com/glennhickey/hal\">HAL tools package</a>.\n") 
    f.write("</P>\n")
    f.write("\n")

def writeLiftoverDocs_credits(f):
    f.write("<H2>Credits</H2>\n")
    f.write("<P>\n")
    f.write("Glenn Hickey, Ngan Nguyen, Joel Armstrong, Benedict Paten.\n")
    f.write("</P>\n")
    f.write("\n")

def writeLiftoverDocs_references(f):
    f.write("<H2>References</H2>\n")
    f.write("<P>\n")
    f.write("Hickey <i>et al.</i>.\n")
    f.write("<a href=\"http://bioinformatics.oxfordjournals.org/content/29/10/1341.short\">HAL: a hierarchical format for storing and analyzing multiple genome alignments.</a>\n")
    f.write("<em>Bioinformatics</em>. 2013 May;29(10):1341-1342.\n")
    f.write("</P>\n")
    f.write("\n")

    f.write("<P>\n")
    f.write("Comparative Assembly Hubs: Web Accessible Browsers for Comparative Genomics\n")
    f.write("</P>\n")
    f.write("\n")

def writeLiftoverDocs(f):
    f.write("<H1>Lifted-over Annotations</H1>\n")
    writeLiftoverDocs_description(f)
    writeLiftoverDocs_methods(f)
    writeLiftoverDocs_credits(f)
    writeLiftoverDocs_references(f)

def makeHubCentralDocs(outdir):
    outfile = os.path.join(outdir, "hubCentral.html")
    f = open(outfile, 'w')
    writeSnakeDocs(f)
    writeLiftoverDocs(f)
    f.close()

def main():
    makeHubCentralDocs(sys.argv[1])

if __name__ == '__main__':
    main()


