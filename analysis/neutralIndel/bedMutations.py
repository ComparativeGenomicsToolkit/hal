#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python3

"""Scan a BED file of mutations
"""
import sys
import os

class BedMutations(object):

    # the official location of these definitions is in
    # hal/mutations/src/halBranchMutations.cpp
    inversionBedTag = "V"
    insertionBedTag = "I"
    deletionBedTag = "D"
    deletionBreakBedTag = "DB"
    transpositionBedTag = "P"
    duplicationBedTag = "U"
    gapInsertionBedTag = "GI"
    gapDeletionBedTag = "GD"
    gapDeletionBreakBedTag = "GDB"
    substitutionBedTag = "S"

    # keep lenient by default, everything but dupes and transpositions
    defaultEvents = [insertionBedTag, gapInsertionBedTag,
                     deletionBedTag, deletionBreakBedTag,
                     gapDeletionBedTag, gapDeletionBreakBedTag]
    
    def __init__(self):
        pass

    # read bed file line by line, storing relevant info in class members
    def scan(self, bedPath, events = None):
        if events is None:
            events = self.defaultEvents
        bedFile = open(bedPath, "r")
        self.sequence = None
        self.prevSequence = None
        self.range = None
        self.prevRange = None
        self.op = None
        self.prevOp = None
        self.events = set(events)
        self.genome = None
        self.ancGenome = None
        
        for line in bedFile:
            tokens = line.split()
            if len(tokens) == 0 or tokens[0][0] == "#":
                continue
            assert len(tokens) >= 6
            if not self.__testIgnore(tokens[3]):
                assert self.genome is None or tokens[5] == self.genome
                self.genome = tokens[5]
                self.ancGenome = tokens[4]
                self.prevSequence = self.sequence
                self.prevRange = self.range
                self.prevOp = self.op
                self.sequence = tokens[0]
                self.range = (int(tokens[1]), int(tokens[2]))
                self.op = tokens[3]
                if self.prevSequence != self.sequence:
                    self.prevSequence = None
                    self.prevRange = None
                    self.prevOp = None
                yield

    # return distance from previous operation
    def distance(self):
        if self.prevRange is None or self.range is None:
            return None
        d = self.range[0] - self.prevRange[1]
        if d < 0:
            raise RuntimeError("Distance between (%d,%d) and (%d,%d) is negative which probably means the mutations bed file is not sorted.  Please regenerate the mutations with halTreeMutations.py (without the --noSort option!)" % (
                self.prevRange[0], self.prevRange[1], self.range[0],
                self.range[1]))
        return d

    def __testIgnore(self, token):
        if token.find(self.substitutionBedTag) == 0:
            return self.substitutionBedTag not in self.events
        return token not in self.events
