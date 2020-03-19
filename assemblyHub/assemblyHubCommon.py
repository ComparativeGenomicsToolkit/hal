#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Objects & functions that are used by multiple assemblyHub modules
"""
import os, re, time
from sonLib.bioio import system  
from toil.job import Job

class MakeAnnotationTracks( Job ):
    def __init__(self, options, outdir, halfile, genome2seq2len, type):
        Job.__init__(self)
        self.options = options
        self.outdir = outdir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.type = type #type has to be bed or wiggle

    def run(self, fileStore):
        #Lift-over Bed/Wiggle annotations of each sample to all other samples
        if self.type == 'bed':
            indirs = self.options.beddirs
        elif self.type == 'bed2':
            indirs = self.options.beddirs2
        else:
            indirs = self.options.wigdirs

        if indirs:
            for indir in indirs:
                annotation = os.path.basename(indir)
                bigdir = os.path.join(self.outdir, "liftover%s" %self.type, annotation)
                if self.type == "wig":
                    from hal.assemblyHub.wigTrack import LiftoverWigFiles
                    self.addChild( LiftoverWigFiles(indir, self.halfile, self.genome2seq2len, bigdir, self.options.noWigLiftover, self.outdir) )
                else:
                    from hal.assemblyHub.bedTrack import LiftoverBedFiles
                    self.addChild( LiftoverBedFiles(indir, self.halfile, self.genome2seq2len, bigdir, self.options.noBedLiftover, self.options.tabbed, self.outdir, self.options) )

class CleanupFiles(Job):
    def __init__(self, files):
        Job.__init__(self)
        self.files = files

    def run(self, fileStore):
        if len(self.files) > 0:
            system( "rm -f %s" % " ".join(self.files) ) #cleanup

def preprocessAnnotationInputs(options, outdir, type):
    bigdirs = []
    if type == 'bed':
        indirs = options.beddirs
    elif type == 'bed2':
        indirs = options.beddirs2
    else:
        indirs = options.wigdirs

    if indirs:
        for indir in indirs:
            annotation = os.path.basename(indir)
            bigdir = os.path.join(outdir, "liftover%s" %type, annotation)
            system("mkdir -p %s" % bigdir)
            bigdirs.append(bigdir)
    #Add the 'final' Bed annotation tracks, no liftover
    bbdirs = options.bbdirs
    if type == 'bed2':
        bbdirs = options.bbdirs2
    elif type == 'wig':
        bbdirs = options.bwdirs
    if bbdirs:
        liftoverdir = os.path.join(outdir, "liftover%s" %type)
        if not os.path.exists(liftoverdir):
            system("mkdir -p %s" %liftoverdir)
        for bbdir in bbdirs:
            annotation = os.path.basename(bbdir)
            bigbeddir = os.path.join(outdir, "liftover%s" %type, annotation)
            bigdirs.append(bigbeddir)
            #Copy bb files to bigbeddir
            if os.path.abspath(bigbeddir) != os.path.abspath(bbdir):
                system("ln -s %s %s" %(os.path.abspath(bbdir), os.path.abspath(bigbeddir)))
                #system("cp -r -L --copy-contents %s %s" %(bbdir, bigbeddir))
    if type == 'bed':
        options.bigbeddirs = bigdirs
    elif type == 'bed2':
        options.bigbeddirs2 = bigdirs
    else:
        options.bigwigdirs = bigdirs

def getFilesByExt(indir, fileExtension):
    #Return all files in indir that have "fileExtension"
    allfiles = os.listdir(indir)
    files = []
    for f in allfiles:
        if os.path.splitext(f) == ".%s" %fileExtension :
            files.append(f)
    return files

def getProperName(name, properName):
    if not name:
        return "NoName"
    newname = name
    if name in properName:
        newname = properName[name]
    items = newname.split()
    return "_".join(items)

def sortByProperName(names, properName):
    proper2names = {}
    for n in names:
        n2 = getProperName(n, properName)
        if n2 not in proper2names:
            proper2names[n2] = [n]
        else:
            proper2names[n2].append(n)
    sortedNames = []
    for n2 in sorted(proper2names.keys()):
        sortedNames.extend( proper2names[n2] )
    return sortedNames



