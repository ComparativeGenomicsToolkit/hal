#!/usr/bin/env python

#Copyright (C) 2013 by Ngan Nguyen
#
#Released under the MIT license, see LICENSE.txt

"""
Objects & functions that are used by multiple assemblyHub modules
"""
import os, re, time
from sonLib.bioio import system  
from jobTree.scriptTree.target import Target

class MakeAnnotationTracks( Target ):
    def __init__(self, options, outdir, halfile, genome2seq2len, type):
        Target.__init__(self)
        self.options = options
        self.outdir = outdir
        self.halfile = halfile
        self.genome2seq2len = genome2seq2len
        self.type = type #type has to be bed or wiggle

    def run(self):
        #Lift-over Bed/Wiggle annotations of each sample to all other samples
        if self.type == 'bed':
            indirs = self.options.beddirs
        else:
            indirs = self.options.wigdirs

        if indirs:
            for indir in indirs:
                annotation = os.path.basename(indir)
                bigdir = os.path.join(self.outdir, "liftover%ss" %self.type, annotation)
                if self.type == "bed":
                    from hal.assemblyHub.bedTrack import LiftoverBedFiles
                    self.addChildTarget( LiftoverBedFiles(indir, self.halfile, self.genome2seq2len, bigdir, self.options.noBedLiftover, self.outdir) )
                else:
                    from hal.assemblyHub.wigTrack import LiftoverWigFiles
                    self.addChildTarget( LiftoverWigFiles(indir, self.halfile, self.genome2seq2len, bigdir, self.options.noWigLiftover, self.outdir) )

class CleanupFiles(Target):
    def __init__(self, files):
        Target.__init__(self)
        self.files = files

    def run(self):
        if len(self.files) > 0:
            system( "rm %s" % " ".join(self.files) ) #cleanup

def preprocessAnnotationInputs(options, outdir, type):
    bigdirs = []
    if type == 'bed':
        indirs = options.beddirs
    else:
        indirs = options.wigdirs

    if indirs:
        for indir in indirs:
            annotation = os.path.basename(indir)
            bigdir = os.path.join(outdir, "liftover%ss" %type, annotation)
            system("mkdir -p %s" % bigdir)
            bigdirs.append(bigdir)
    #Add the 'final' Bed annotation tracks, no liftover
    bbdirs = options.bbdirs
    if type == 'wig':
        bbdirs = options.bwdirs
    if bbdirs:
        liftoverdir = os.path.join(outdir, "liftover%ss" %type)
        if not os.path.exists(liftoverdir):
            system("mkdir -p %s" %liftoverdir)
        for bbdir in bbdirs:
            annotation = os.path.basename(bbdir)
            bigbeddir = os.path.join(outdir, "liftover%ss" %type, annotation)
            bigdirs.append(bigbeddir)
            #Copy bb files to bigbeddir
            if os.path.abspath(bigbeddir) != os.path.abspath(bbdir):
                system("cp -r -L --copy-contents %s %s" %(bbdir, bigbeddir))
    if type == 'bed':
        options.bigbeddirs = bigdirs
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
    newname = name
    if name in properName:
        newname = properName[name]
    items = newname.split()
    return "_".join(items)


