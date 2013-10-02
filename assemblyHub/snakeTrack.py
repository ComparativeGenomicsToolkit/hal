#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""Snake tracks
"""
import re

def writeTrackDb_snakes(f, halfile, genomes, currgenome, properName):
    for i, genome in enumerate(genomes):
        if re.search(genome, currgenome): #current genome
            continue
        #SNAKE TRACKS
        genomeProperName = genome
        if genome in properName:
            genomeProperName = properName[genome]
        f.write("track snake%s\n" %genome)
        f.write("longLabel %s\n" %genomeProperName)
        f.write("shortLabel %s\n" %genomeProperName)
        f.write("otherSpecies %s\n" %genome)
        f.write("visibility full\n")
        f.write("priority %d\n" %(i + 2))
        f.write("bigDataUrl %s\n" % halfile)
        f.write("type halSnake\n")
        f.write("group snake\n")
        f.write("\n")

