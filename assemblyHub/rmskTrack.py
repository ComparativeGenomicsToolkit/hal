#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""RepeatMasker track
"""
import os
from sonLib.bioio import system  
from optparse import OptionGroup

def writeTrackDb_rmsk(f, rmskdir, genomedir):
    if not os.path.exists(rmskdir):
        return
    f.write("track repeatMasker\n")
    f.write("compositeTrack on\n")
    f.write("shortLabel RepeatMasker\n")
    f.write("longLabel Repeating Elements by RepeatMasker\n")
    f.write("group map\n")
    f.write("visibility dense\n")
    f.write("type bed 3 .\n")
    f.write("noInherit on\n")
    f.write("html ../documentation/repeatMasker\n")
    f.write("\n")
    
    system("ln -s %s %s" %(os.path.abspath(rmskdir), os.path.join(genomedir, "repeatMasker")))
    files = os.listdir(rmskdir)
    for i, file in enumerate(files):
        element = file.split('.')[0]
        f.write("\ttrack repeatMasker%s\n" %element)
        f.write("\tparent repeatMasker\n")
        f.write("\tshortLabel %s\n" %element)
        f.write("\tlongLabel %s Repeating Elements by RepeatMasker\n" %element)
        f.write("\tpriority %d\n" %i)
        f.write("\tspectrum on\n")
        f.write("\tmaxWindowToDraw 10000000\n")
        f.write("\tcolorByStrand 50,50,150 150,50,50\n")
        f.write("\ttype bigBed 6 +\n")
        f.write("\tbigDataUrl repeatMasker/%s.bb\n" %element)
        f.write("\n")
    f.write("\n")

def addRmskOptions(parser):
    #group = parser.add_argument_group("REPEATMASKER", "RepeatMasker options.")
    group = parser.add_argument_group("REPEATMASKER")
    group.add_argument('--rmskDir', dest='rmskdir', help="Directory containing repeatMasker's output files for each genome. Format: rmskDir/ then genome1/ then genome.rmsk.SINE.bb, genome.rmsk.LINE.bb, ... ")
    group = parser.add_argument_group(group)

def checkRmskOptions(parser, options):
    if options.rmskdir:
        if not os.path.exists(options.rmskdir) or not os.path.isdir(options.rmskdir):
            parser.error("RepeatMasker directory %s does not exist or is not a directory.\n" %options.rmskdir)





