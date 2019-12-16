#!/usr/bin/env python3
"""Runs ancestorsML on all unique columns that contain at least one
ancestor.
"""
import math
from argparse import ArgumentParser
from collections import defaultdict

from sonLib.bioio import getTempFile, system, popenCatch
from sonLib.nxnewick import NXNewick
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target

class Setup(Target):
    def __init__(self, halFile, phyloPModel, jobsPerGenome, threshold):
        Target.__init__(self)
        self.halFile = halFile
        self.phyloPModel = phyloPModel
        self.jobsPerGenome = jobsPerGenome
        self.threshold = threshold

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
        self.setFollowOnTarget(RunAncestorsMLParallel(self.halFile, self.phyloPModel, bedFiles, self.jobsPerGenome, self.threshold))

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
    def __init__(self, halFile, phyloPModel, bedFileDict, jobsPerGenome, threshold):
        Target.__init__(self)
        self.halFile = halFile
        self.phyloPModel = phyloPModel
        self.bedFileDict = bedFileDict
        self.jobsPerGenome = jobsPerGenome
        self.threshold = threshold

    def run(self):
        outputsPerGenome = {}
        for genome, bedFile in list(self.bedFileDict.items()):
            outputsPerGenome[genome] = []
            with open(bedFile) as f:
                numLines = sum(1 for line in f)
            linesPerJob = int(math.ceil(float(numLines)/self.jobsPerGenome))
            if linesPerJob == 0:
                linesPerJob = 1
            for start in range(0, numLines, linesPerJob):
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
                                                   self.threshold, output))
                outputsPerGenome[genome].append(output)
        self.setFollowOnTarget(WriteNucleotides(outputsPerGenome, self.halFile))

class RunAncestorsML(Target):
    def __init__(self, halFile, genome, bedForJob, modelFile, threshold, output):
        Target.__init__(self)
        self.halFile = halFile
        self.genome = genome
        self.bedForJob = bedForJob
        self.modelFile = modelFile
        self.threshold = threshold
        self.output = output

    def run(self):
        system("ancestorsML --printWrites --bed %s --thresholdN %s %s %s %s > %s" % (self.bedForJob, self.threshold, self.halFile, self.genome, self.modelFile, self.output))

class WriteNucleotides(Target):
    def __init__(self, inputsPerGenome, halFile):
        Target.__init__(self)
        self.inputsPerGenome = inputsPerGenome
        self.halFile = halFile

    def run(self):
        counts = defaultdict(int)
        nCounts = defaultdict(int)
        for genome, inputs in list(self.inputsPerGenome.items()):
            for input in inputs:
                with open(input) as f:
                    for line in f:
                        genome, position, old, new = line.strip().split('\t')
                        counts[genome] += 1
                        if new == 'N':
                            nCounts[genome] += 1

                system("halWriteNucleotides %s %s" % (self.halFile, input))

        for genome in counts:
            self.logToMaster('Changed %s nucleotides of genome %s (%s to N)' % (counts[genome], genome, nCounts[genome]))

if __name__ == '__main__':
    from ancestorsMLMP import * # required for jobTree
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('halPath', help='hal file (will be modified, make a backup!)')
    parser.add_argument('phyloPModel', help='phyloP model file generated using halPhyloPTrain.py')
    parser.add_argument('--jobsPerGenome', help='maximum number of jobs per genome',
                        type=int, default=2000)
    parser.add_argument('--threshold', help='bases with a posterior probability less '
                        'than this will be set to N (set to 0.0 to disable)', type=float,
                        default=0.9)
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()
    Stack(Setup(opts.halPath, opts.phyloPModel, opts.jobsPerGenome, opts.threshold)).startJobTree(opts)
