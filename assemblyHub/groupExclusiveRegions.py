#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Creating bed tracks of group exclusive regions
"""
import os
from sonLib.bioio import system  
from toil.job import Job
from optparse import OptionGroup
from hal.assemblyHub.treeCommon import iterAllClades 

#---------- BRANCH-SPECIFIC REGIONS (regions shared by a branch and not by anybody else in the tree ----------
class GetCladeExclusiveRegions(Job):
    #outdir: liftoverbeds/CladeExclusive
    def __init__(self, halfile, tree, bigbeddir, maxOut, minIn):
        Job.__init__(self)
        self.halfile = halfile
        self.tree = tree
        self.bigbeddir = bigbeddir
        self.maxOut = maxOut
        self.minIn = minIn

    def run(self, fileStore):
        allClades = iterAllClades(self.tree.root)
        outdir = os.path.join(self.bigbeddir, "CladeExclusive")
        system("mkdir -p %s" %(outdir))
        
        for clade in allClades:
            if len(clade) == 0:
                continue
            self.addChild( GetCladeExclusive(self.halfile, clade, outdir, self.maxOut, self.minIn) )

class GetCladeExclusive(Job):
    def __init__(self, halfile, names, outdir, maxOut, minIn):
        Job.__init__(self)
        self.halfile = halfile
        self.names = names
        self.outdir = outdir
        self.maxOut = maxOut
        self.minIn = minIn

    def run(self, fileStore):
        cladedir = os.path.join(self.outdir, self.names[0])
        print(self.names)
        print(cladedir)
        system("mkdir -p %s" %cladedir)
        minIn = self.minIn
        if not minIn or minIn > len(self.names):
            minIn = len(self.names)
        maxOut = self.maxOut
        if maxOut > minIn:
            maxOut = minIn
        
        #Get clade exclusive regions, w.r.t root
        outbed = os.path.join(cladedir, "%s.bed" %self.names[0])
        cmd = "findRegionsExclusivelyInGroup --maxOutgroupGenomes %d --minIngroupGenomes %d %s %s %s > %s" \
              %(maxOut, minIn, self.halfile, self.names[0], ",".join(self.names), outbed )
        system(cmd)
        
        #Convert to bigbed
        chrsizefile = os.path.join(self.outdir, "../..", self.names[0], "chrom.sizes")
        outbigbed = os.path.join(cladedir, "%s.bb" %self.names[0])
        system("bedToBigBed %s %s %s" %(outbed, chrsizefile, outbigbed))

        #Lift over to the children genomes:
        for child in self.names[1:]:
            self.addChild(LiftoverCladeExclusive(cladedir, self.halfile, self.names[0], outbed, child, chrsizefile))
        self.addFollowOn(CleanupCladeExclusive(cladedir))

class CleanupCladeExclusive(Job):
    def __init__(self, cladedir):
        Job.__init__(self)
        self.cladedir = cladedir

    def run(self, fileStore):
        system("rm %s/*bed" %self.cladedir)

class LiftoverCladeExclusive(Job):
    def __init_(self, cladedir, halfile, query, queryBed, target, chrsizefile):
        Job.__init__(self)
        self.cladedir = cladedir
        self.halfile = halfile
        self.query = query
        self.queryBed = queryBed
        self.target = target
        self.chrsizefile = chrsizefile

    def run(self, fileStore):
        bedfile = os.path.join(self.cladedir, "%s.bed" %self.target)
        system("halLiftover %s %s %s %s %s" %(self.halfile, self.query, self.queryBed, self.target, bedfile))
        #Convert to big bed:
        bigbedfile = os.path.join(self.cladedir, "%s.bb" %self.target)
        system("bedToBigBed %s %s %s" %(bedfile, self.chrsizefile, bigbedfile))

def addExclusiveRegionOptions(parser):
    group = parser.add_argument_group("CLADE EXCLUSIVE REGIONS", "Requirements of regions that are exclusive to subgroup of genomes.")
    group.add_argument('--cladeExclusiveRegions', dest='cladeExclusive', action='store_true', default=False, help='If specified, will generate tracks of regions that are exclusive to each branch (including leaf "branches", which will be genome-exclusive regions) on the tree. ')
    group.add_argument('--maxOutgroupGenomes', dest='maxOut', type=int, default=0, help='Maximum number of outgroup genomes that a region is allowed to be in. ')
    group.add_argument('--minIngroupGenomes', dest='minIn', type=int, help='Minimum number of ingroup genomes that a region must appear in. Default=all ingroup genomes (branch node and all its children).')
    group = parser.add_argument_group(group)







