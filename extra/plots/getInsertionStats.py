#!/usr/bin/env python
from argparse import ArgumentParser
from sonLib.bioio import popenCatch, system, getTempFile
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target

class Setup(Target):
    def __init__(self, opts):
        Target.__init__(self)
        self.opts = opts

    def run(self):
        genomes = popenCatch("halStats --genomes %s" % self.opts.halPath).split()
        outputs = []
        for genome in genomes:
            tempOutput = getTempFile(rootDir=self.getGlobalTempDir())
            outputs.append(tempOutput)
            self.addChildTarget(ExtractInsertions(self.opts.halPath, genome, tempOutput, self.opts.samplePerGenome, self.opts.samples, self.opts.noGaps))
        self.setFollowOnTarget(ReduceOutputs(outputs, self.opts.output, not self.opts.samplePerGenome, self.opts.samples))

class ExtractInsertions(Target):
    def __init__(self, halPath, genome, output, doSampling, numSamples, removeGaps):
        Target.__init__(self)
        self.halPath = halPath
        self.genome = genome
        self.output = output
        self.doSampling = doSampling
        self.numSamples = numSamples
        self.removeGaps = removeGaps

    def run(self):
        insertionBed = getTempFile(rootDir=self.getGlobalTempDir())
        system("halAlignedExtract --complement %s %s > %s" % (self.halPath, self.genome, insertionBed))
        if self.doSampling:
            # Sample per-genome instead of overall
            temp = getTempFile(rootDir=self.getGlobalTempDir())
            system("shuf %s | head -n %d > %s" % (insertionBed, self.numSamples, temp))
            system("mv %s %s" % (temp, insertionBed))
        outFile = open(self.output, 'w')
        for line in open(insertionBed):
            fields = line.split()
            if len(fields) >= 3:
                seq = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if self.removeGaps:
                    # Get rid of gaps. This will make everything slow
                    # as hell and we should really load the sequence
                    # once and then query it. but whatever
                    dna = popenCatch("hal2fasta %s %s --sequence %s --start %d --length %d" % (self.halPath, seq, start, end - start))
                    # remove header
                    dna = dna.split("\n")[1:]
                    if 'N' in dna or 'n' in dna:
                        # Found a gap
                        continue
                outFile.write("%d\t%s\t%s\n" % (end - start, self.genome, seq))
        outFile.close()

class ReduceOutputs(Target):
    def __init__(self, outputs, outPath, doSampling, numSamples):
        Target.__init__(self)
        self.outputs = outputs
        self.outPath = outPath
        self.doSampling = doSampling
        self.numSamples = numSamples

    def run(self):
        outFile = open(self.outPath, 'w')
        for output in self.outputs:
            for line in open(output):
                outFile.write(line)
        if self.doSampling:
            # Sample overall instead of per-genome
            temp = getTempFile(rootDir=self.getGlobalTempDir())
            system("shuf %s | head -n %d > %s" % (self.outPath, self.numSamples, temp))
            system("mv %s %s" % (temp, self.outPath))

if __name__ == '__main__':
    from getInsertionStats import * # required for jobTree
    parser = ArgumentParser()
    parser.add_argument('halPath', help='hal file')
    parser.add_argument('output', help='output tsv')
    parser.add_argument('--samples', help='number of samples per genome', default=10000, type=int)
    parser.add_argument('--samplePerGenome', action="store_true", help="sample n samples per genome instead of n samples overall")
    parser.add_argument('--noGaps', action="store_true", help="remove any sequences with Ns in them")
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()
    Stack(Setup(opts)).startJobTree(opts)
