#!/usr/bin/env python
"""Runs ancestorsML on all unique columns that contain at least one
ancestor.
"""
import math
from argparse import ArgumentParser
from sonLib.bioio import getTempFile, system, popenCatch
from sonLib.nxnewick import NXNewick
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target

class Setup(Target):
    def __init__(self, halFile, phyloPModel, jobsPerGenome):
        Target.__init__(self)
        self.halFile = halFile
        self.phyloPModel = phyloPModel
        self.jobsPerGenome = jobsPerGenome

    def run(self):
        # Find all ancestral genomes using the tree.
        newickStr = popenCatch("halStats --tree %s" % self.halFile)
        tree = NXNewick().parseString(newickStr)
        bedFiles = {} # genome => bed files of inserted columns
        for nodeId in tree.postOrderTraversal():
            if len(tree.getChildren(nodeId)) == 0:
                # leaf node, skip
                continue
            assert tree.hasName(nodeId)
            genome = tree.getName(nodeId)
            bedFileForGenome = getTempFile(rootDir=self.getGlobalTempDir())
            bedFiles[genome] = bedFileForGenome
            self.addChildTarget(GetInsertedColumnBed(self.halFile, genome, bedFileForGenome))
        self.setFollowOnTarget(RunAncestorsMLParallel(self.halFile, self.phyloPModel, bedFiles, self.jobsPerGenome))

class GetInsertedColumnBed(Target):
    """Gets a bed file containing all columns inserted in this genome."""
    def __init__(self, halFile, genome, outputPath):
        Target.__init__(self)
        self.halFile = halFile
        self.genome = genome
        self.outputPath = outputPath

    def run(self):
        system("halAlignedExtract --complement %s %s > %s" % (self.halFile,
                                                              self.genome,
                                                              self.outputPath))

class RunAncestorsMLParallel(Target):
    def __init__(self, halFile, phyloPModel, bedFileDict, jobsPerGenome):
        Target.__init__(self)
        self.halFile = halFile
        self.phyloPModel = phyloPModel
        self.bedFileDict = bedFileDict
        self.jobsPerGenome = jobsPerGenome

    def run(self):
        outputsPerGenome = {}
        for genome, bedFile in self.bedFileDict.items():
            outputsPerGenome[genome] = []
            numLines = int(popenCatch("wc -l %s | cut -d' ' -f 1" % bedFile))
            linesPerJob = int(math.ceil(float(numLines)/self.jobsPerGenome))
            if linesPerJob == 0:
                linesPerJob = 1
            for start in xrange(0, numLines, linesPerJob):
                end = start + linesPerJob
                if end > numLines:
                    end = numLines
                bedForJob = getTempFile(rootDir=self.getGlobalTempDir())
                system("head -n %d %s | tail -n %d > %s" % (start + linesPerJob,
                                                            bedFile,
                                                            end - start,
                                                            bedForJob))
                output = getTempFile(rootDir=self.getGlobalTempDir())
                self.addChildTarget(RunAncestorsML(self.halFile, genome,
                                                   bedForJob, self.phyloPModel,
                                                   output))
                outputsPerGenome[genome].append(output)
        self.setFollowOnTarget(WriteNucleotides(outputsPerGenome, self.halFile))

class RunAncestorsML(Target):
    def __init__(self, halFile, genome, bedForJob, modelFile, output):
        Target.__init__(self)
        self.halFile = halFile
        self.genome = genome
        self.bedForJob = bedForJob
        self.modelFile = modelFile
        self.output = output

    def run(self):
        system("ancestorsML --printWrites --bed %s %s %s %s > %s" % (self.bedForJob, self.halFile, self.genome, self.modelFile, self.output))

class WriteNucleotides(Target):
    def __init__(self, inputsPerGenome, halFile):
        Target.__init__(self)
        self.inputsPerGenome = inputsPerGenome
        self.halFile = halFile

    def run(self):
        for genome, inputs in self.inputsPerGenome.items():
            for input in inputs:
                system("halWriteNucleotides %s %s" % (self.halFile, input))

if __name__ == '__main__':
    from ancestorsMLMP import * # required for jobTree
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('halPath', help='hal file (will be modified, make a backup!)')
    parser.add_argument('phyloPModel', help='phyloP model file generated using halPhyloPTrain.py')
    parser.add_argument('--jobsPerGenome', help='maximum number of jobs per genome',
                        type=int, default=2000)
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()
    Stack(Setup(opts.halPath, opts.phyloPModel, opts.jobsPerGenome)).startJobTree(opts)
