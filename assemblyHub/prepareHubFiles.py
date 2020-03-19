#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Make "hub.txt", "groups.txt", files that are required by AssemblyHub
Also prepare description.html files   
"""

import os, sys
from sonLib.bioio import system  
from optparse import OptionGroup
from hal.assemblyHub.assemblyHubCommon import getProperName
from Bio import Phylo
from hal.assemblyHub.treeCommon import isBinaryTree

def writeDescriptionFile(genome, outdir):
    filename = os.path.join(outdir, "description.html")
    f = open(filename, 'w')
    f.write("%s\n" %genome)
    f.close()
    return

def writeTrackDb_composite_html(file, treeFile):
    f = open(file, 'w')
    #HACK:
    #huburl = "http://hgwdev.cse.ucsc.edu/~nknguyen/ecoli/hub/TEST2"
    huburl = "http://hgwdev.cse.ucsc.edu/~nknguyen/birds/birds2"
    basename = os.path.basename(treeFile)
    f.write("<img src=\"%s/%s\">\n" %(huburl, basename))
    f.close()

def writeTrackDb_compositeStart(f, shortLabel, longLabel, bbdirs, bwdirs, genomes, properName, url, img):
    #Composite track includes all annotations in BED & WIGGLE formats, their lifted-over tracks, and Snake tracks
    f.write("track hubCentral\n")
    f.write("compositeTrack on\n")
    f.write("shortLabel %s\n" %shortLabel)
    f.write("longLabel %s\n" %longLabel)
    f.write("group comphub\n")

    bedtracktypes = [os.path.basename(b.rstrip('/')) for b in bbdirs]
    bedstr = " ".join(["%s=%s" %(item, item) for item in bedtracktypes])
    wigtracktypes = [os.path.basename(b.rstrip('/')) for b in bwdirs]
    wigstr = " ".join(["%s=%s" %(item, item) for item in wigtracktypes])
    
    f.write("subGroup1 view Track_Type Snake=Alignments %s %s\n" %(bedstr, wigstr))
    
    genomeStr = " ".join(["%s=%s" %(g, getProperName(g, properName)) for g in genomes])
    f.write("subGroup2 orgs Organisms %s\n" %genomeStr)
    f.write("dragAndDrop subTracks\n")
    f.write("#allButtonPair on\n")
    #f.write("sortOrder view=+ orgs=+\n")
    f.write("dimensions dimensionX=view dimensionY=orgs\n") 
    f.write("noInherit on\n")
    f.write("priority 0\n")
    f.write("centerLabelsDense on\n")
    f.write("visibility full\n")
    f.write("html ../documentation/hubCentral\n")

    if url and img:
        imgurl = os.path.join(url, os.path.basename(img))
        f.write("treeImage %s\n" %imgurl)

    f.write("type bigBed 3\n")
    f.write("\n")

def writeTrackDb_compositeSubTrack(f, name, visibility):
    f.write("\ttrack hubCentral%s\n" %name)
    f.write("\tshortLabel %s\n" %name) 
    f.write("\tview %s\n" %name)
    f.write("\tvisibility %s\n" %visibility)
    f.write("\tsubTrack hubCentral\n")
    f.write("\n")

def writeGroupFile(outdir, hubLabel, annotations):
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

    f.write("name comphub\n")
    f.write("label %s\n" % hubLabel)
    f.write("priority 3\n")
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

#=========== READ FILES ===========
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

#=========== OPTIONS =============
def addHubOptions(parser):
    group = parser.add_argument_group("HUB INFORMATION")
    group.add_argument('--hub', dest='hubLabel', default='myHub', help='a single-word name of the directory containing the track hub files. Not displayed to hub users. ')
    group.add_argument('--shortLabel', dest='shortLabel', default='my hub', help='the short name for the track hub. Suggested maximum length is 17 characters. Displayed as the hub name on the Track Hubs page and the track group name on the browser tracks page. ')
    group.add_argument('--longLabel', dest='longLabel', default='my hub', help='a longer descriptive label for the track hub. Suggested maximum length is 80 characters. Displayed in the description field on the Track Hubs page. ')
    group.add_argument('--email', dest='email', default='NoEmail', help='the contact to whom questions regarding the track hub should be directed. ')
    group.add_argument('--genomes', dest='genomes', help='File specified list of genomes to make browser for. If specified, only create browsers for these genomes in the order provided by the list. Otherwise create browsers for all genomes in the input hal file')
    group.add_argument('--rename', dest='rename', help='File that maps halfile genomeNames to names displayed on the browser. Format: <halGenomeName>\\t<genomeNameToDisplayOnBrowser>. ') 
    group.add_argument('--tree', dest='treeFile', help='Newick binary tree. The order of the tracks and the default track layout will be based on this tree if option "genomes" is not specified. If not specified, try to extract the newick tree from the input halfile.')
    group.add_argument('--url', dest='url', help='Public url of the hub location')
    group.add_argument('--twobitdir', dest='twobitdir', help='Optional. Directory containing the 2bit files of each genomes. Default: extract from the input hal file.')
    group = parser.add_argument_group(group)

def checkHubOptions(parser, options):
    if options.genomes:
        options.genomes = readList(options.genomes)
    options.properName = {}
    if options.rename and os.path.exists(options.rename):
        options.properName = readRename(options.rename)
    
    options.treeFig = None
    options.leaves = None
    options.tree = None
    if options.treeFile and not os.path.exists(options.treeFile):
        parser.error("The tree file %s does not exist.\n" %options.tree)
    elif options.treeFile:
        tree = Phylo.read(options.treeFile, 'newick')
        if isBinaryTree(tree):
            options.tree = tree
        else:
            sys.stderr.write("Warnning: tree %s is not a binary tree. Will be ignored!" %options.treeFile)

