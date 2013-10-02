#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""
Make "hub.txt", "groups.txt", files that are required by AssemblyHub
Also prepare description.html files   
"""

import os
from sonLib.bioio import system  
from optparse import OptionGroup

def writeDescriptionFile(genome, outdir):
    filename = os.path.join(outdir, "description.html")
    f = open(filename, 'w')
    f.write("%s\n" %genome)
    f.close()
    return

def writeGroupFile(outdir, annotations):
    filename = os.path.join(outdir, "groups.txt")
    f = open(filename, 'w')
    f.write("name user\n")
    f.write("label Custom\n")
    f.write("priority 1\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")

    f.write("name map\n")
    f.write("label Mapping\n")
    f.write("priority 2\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")

    f.write("name snake\n")
    f.write("label Alignment Snakes\n")
    f.write("priority 3\n")
    f.write("defaultIsClosed 0\n")
    f.write("\n")
    
    for annotation in annotations:
        f.write("name annotation%s\n" %annotation)
        f.write("label %s Annotations\n" % annotation.capitalize() )
        f.write("priority 3\n")
        f.write("defaultIsClosed 1\n")
        f.write("\n")   

    f.write("name exp\n")
    f.write("label Experimental\n")
    f.write("priority 4\n")
    f.write("defaultIsClosed 1\n")
    f.write("\n")
    
    f.close()

def writeHubFile(outdir, options):
    hubfile = os.path.join(outdir, "hub.txt")
    f = open(hubfile, "w")
    f.write("hub %s\n" %options.hubLabel)
    f.write("shortLabel %s\n" %options.shortLabel)
    f.write("longLabel %s\n" %options.longLabel)
    f.write("genomesFile genomes.txt\n")
    f.write("email %s\n" %options.email)
    f.close()

def readList(file):
    items = []
    f = open(file, 'r')
    for line in f:
        items.append(line.strip())
    f.close()
    return items

def readRename(file):
    name2new = {}
    f = open(file, 'r')
    for line in f:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        items = line.split('\t')
        if len(items) >=2:
            name2new[items[0]] = items[1]
    f.close()
    return name2new

def addHubOptions(parser):
    group = OptionGroup(parser, "HUB INFORMATION")
    group.add_option('--hub', dest='hubLabel', default='myHub', help='a single-word name of the directory containing the track hub files. Not displayed to hub users. Default=%default')
    group.add_option('--shortLabel', dest='shortLabel', default='my hub', help='the short name for the track hub. Suggested maximum length is 17 characters. Displayed as the hub name on the Track Hubs page and the track group name on the browser tracks page. Default=%default')
    group.add_option('--longLabel', dest='longLabel', default='my hub', help='a longer descriptive label for the track hub. Suggested maximum length is 80 characters. Displayed in the description field on the Track Hubs page. Default=%default')
    group.add_option('--email', dest='email', default='NoEmail', help='the contact to whom questions regarding the track hub should be directed. Default=%default')
    group.add_option('--genomes', dest='genomes', help='File specified list of genomes to make browser for. If specified, only create browsers for these genomes in the order provided by the list. Otherwise create browsers for all genomes in the input hal file')
    group.add_option('--rename', dest='rename', help='File that maps halfile genomeNames to names displayed on the browser. Format: <halGenomeName>\\t<genomeNameToDisplayOnBrowser>. Default=%default') 
    parser.add_option_group(group)

def checkHubOptions(parser, options):
    if options.genomes:
        options.genomes = readList(options.genomes)
    options.properName = {}
    if options.rename and os.path.exists(options.rename):
        options.properName = readRename(options.rename)

