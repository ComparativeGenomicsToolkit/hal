#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""Snake tracks
"""
from optparse import OptionGroup
import re

def addSnakeOptions(parser):
    group = OptionGroup(parser, "SNAKE TRACKS", "Snake track options")
    group.add_option('--selfAlignmentSnakes', dest="selfAlignmentTrack",
                     help="Produce a self-alignment snake track for every genome",
                     action="store_true", default=False)
    parser.add_option_group(group)

def writeTrackDb_snakes(f, halfile, genomes, subgenomes, currgenome, properName, snpwidth=None, doSelfAlignment=False):
    for i, genome in enumerate(genomes):
        if not doSelfAlignment and genome == currgenome: #current genome
            continue
        #SNAKE TRACKS
        genomeProperName = genome
        if genome in properName:
            genomeProperName = properName[genome]
        if genome == currgenome:
            genomeProperName += " (self)"
        f.write("\t\ttrack snake%s\n" %genome)
        f.write("\t\tlongLabel %s\n" %genomeProperName)
        f.write("\t\tshortLabel %s\n" %genomeProperName)
        f.write("\t\totherSpecies %s\n" %genome)
        if genome in subgenomes:
            f.write("\t\tvisibility full\n")
            f.write("\t\tparent hubCentralAlignments\n")
        else:
            f.write("\t\tvisibility hide\n")
            f.write("\t\tparent hubCentralAlignments off\n")
        if snpwidth:
            f.write("\t\tshowSnpWidth %d\n" % snpwidth)
        f.write("\t\tpriority %d\n" %(i + 2))
        f.write("\t\tbigDataUrl %s\n" % halfile)
        f.write("\t\ttype halSnake\n")
        f.write("\t\tgroup snake\n")
        f.write("\t\tsubGroups view=Snake orgs=%s\n" %genome)
        f.write("\n")

