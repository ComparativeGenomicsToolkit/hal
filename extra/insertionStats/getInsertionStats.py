#!/usr/bin/env python3
"""Get a TSV of (insertion size, length, seq, genome) sampled from a HAL
file."""
from argparse import ArgumentParser
from sonLib.bioio import popenCatch, system, getTempFile, fastaRead
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target

def countMaskedBases(string):
    ret = 0
    for char in string:
        if char != 'A' and char != 'C' and char != 'T' and char != 'G':
            ret += 1
    return ret

class Setup(Target):
    def __init__(self, opts):
        Target.__init__(self)
        self.opts = opts

    def run(self):
        genomes = popenCatch("halStats --genomes %s" % self.opts.halPath).split()
        # main outputs, entirely inserted sequence outputs, total inserted bases outputs
        outputss = [[], [], []]
        for genome in genomes:
            # Get a temp file to hold the genome's output, which will
            # be concatenated with the others at the end
            tempOutput = getTempFile(rootDir=self.getGlobalTempDir())
            outputss[0].append(tempOutput)

            # Create a temp file to hold entirely inserted seqs, if needed
            tempEntirelyInsertedSequencesPath = None
            if self.opts.entirelyInsertedSequencesPath is not None:
                tempEntirelyInsertedSequencesPath = getTempFile(rootDir=self.getGlobalTempDir())
                outputss[1].append(tempEntirelyInsertedSequencesPath)

            # Create a temp file to hold total inserted bases, if needed
            tempTotalInsertedBasesPath = None
            if self.opts.totalInsertedBasesPath is not None:
                tempTotalInsertedBasesPath = getTempFile(rootDir=self.getGlobalTempDir())
                outputss[2].append(tempTotalInsertedBasesPath)

            self.addChildTarget(ExtractInsertions(self.opts.halPath, genome, tempOutput, self.opts.samplePerGenome, self.opts.samples, self.opts.noGaps, tempEntirelyInsertedSequencesPath, tempTotalInsertedBasesPath))
        self.setFollowOnTarget(ReduceOutputs(outputss, [self.opts.output, self.opts.entirelyInsertedSequencesPath, self.opts.totalInsertedBasesPath], [not self.opts.samplePerGenome, False, False], [self.opts.samples, None, None], ['insertionSize\tgenome\tseq\tmaskedBases', 'insertionSize\tgenome\tseq\tmaskedBases', 'genome\ttotalInsertedBases']))

class ExtractInsertions(Target):
    def __init__(self, halPath, genome, output, doSampling, numSamples, removeGaps, entirelyInsertedSequencePath, totalInsertedBasesPath):
        Target.__init__(self)
        self.halPath = halPath
        self.genome = genome
        self.output = output
        self.doSampling = doSampling
        self.numSamples = numSamples
        self.removeGaps = removeGaps
        self.entirelyInsertedSequencePath = entirelyInsertedSequencePath
        self.totalInsertedBasesPath = totalInsertedBasesPath

    def getFastaDict(self):
        temp = getTempFile(rootDir=self.getGlobalTempDir())
        system("hal2fasta %s %s > %s" % (self.halPath, self.genome, temp))
        ret = {}
        for header, seq in fastaRead(temp):
            ret[header] = seq
        return ret

    def logEntirelyInsertedSequences(self, fastaDict, chromSizes, insertionBed):
        outFile = open(self.entirelyInsertedSequencePath, 'w')
        for line in open(insertionBed):
            fields = line.split()
            if len(fields) >= 3:
                seq = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if end - start == chromSizes[seq]:
                    dna = fastaDict[seq][start:end]
                    maskedBases = countMaskedBases(dna)
                    outFile.write("%d\t%s\t%s\t%s\n" % (end - start, self.genome, seq, maskedBases))

    def logTotalInsertedBases(self, insertionBed):
        outFile = open(self.totalInsertedBasesPath, 'w')
        total = 0
        for line in open(insertionBed):
            fields = line.split()
            if len(fields) >= 3:
                total += int(fields[2]) - int(fields[1])
        outFile.write('%s\t%d\n' % (self.genome, total))

    def run(self):
        fastaDict = self.getFastaDict()
        chromSizes = dict([(x[0], len(x[1])) for x in list(fastaDict.items())])

        insertionBed = getTempFile(rootDir=self.getGlobalTempDir())
        system("halAlignedExtract --complement %s %s > %s" % (self.halPath, self.genome, insertionBed))

        if self.entirelyInsertedSequencePath is not None:
            # Look for insertions that cover an entire sequence
            self.logEntirelyInsertedSequences(fastaDict, chromSizes, insertionBed)

        if self.totalInsertedBasesPath is not None:
            # Output the total number of inserted bases in this genome
            self.logTotalInsertedBases(insertionBed)

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
                dna = fastaDict[seq][start:end]
                if self.removeGaps:
                    # Get rid of gaps.
                    if 'N' in dna or 'n' in dna:
                        # Found a gap
                        continue
                maskedBases = countMaskedBases(dna)
                outFile.write("%d\t%s\t%s\t%s\n" % (end - start, self.genome, seq, maskedBases))
        outFile.close()

class ReduceOutputs(Target):
    def __init__(self, outputss, outPaths, doSamplings, numSampless, headers):
        Target.__init__(self)
        self.outputss = outputss
        self.outPaths = outPaths
        self.doSamplings = doSamplings
        self.numSampless = numSampless
        self.headers = headers

    def run(self):
        for outputs, outPath, doSampling, numSamples, header in zip(self.outputss, self.outPaths, self.doSamplings, self.numSampless, self.headers):
            if outPath is None:
                # No output for this file
                continue
            outFile = open(outPath, 'w')
            if header is not None:
                outFile.write('%s\n' % (header))
            for output in outputs:
                for line in open(output):
                    outFile.write(line)
            if doSampling:
                # Sample overall instead of per-genome
                temp = getTempFile(rootDir=self.getGlobalTempDir())
                system("shuf %s | head -n %d > %s" % (outPath, numSamples, temp))
                system("mv %s %s" % (temp, outPath))

if __name__ == '__main__':
    from getInsertionStats import * # required for jobTree
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('halPath', help='hal file')
    parser.add_argument('output', help='output tsv for a sample of insertions')
    parser.add_argument('--samples', help='number of samples', default=10000, type=int)
    parser.add_argument('--samplePerGenome', action="store_true", help="sample n samples per genome instead of n samples overall")
    parser.add_argument('--noGaps', action="store_true", help="remove any sequences with Ns in them")
    parser.add_argument('--entirelyInsertedSequencesPath', help="tsv to store information about any sequences that are completely unaligned")
    parser.add_argument('--totalInsertedBasesPath', help="tsv to store total inserted bases per-genome")
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()
    Stack(Setup(opts)).startJobTree(opts)
