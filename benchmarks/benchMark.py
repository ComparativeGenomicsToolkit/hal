#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt


import argparse
import os
import sys
import traceback
import time

from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import popenCatch
from sonLib.bioio import system


def runExecutable(command):
    cmdLine = "./runAndGetResources.py \'%s\'" % command
    output = popenCatch(cmdLine)
    values = output[1:-2].split(',')
    return (float(values[0]), float(values[1]), float(values[2]))
    
def runHalGen(preset, hdf5Chunk, hdf5Compression, outPath):
    return runExecutable("halRandGen --preset %s --hdf5Chunk %d\
    --hdf5Compression %d %s" % (preset, hdf5Chunk, hdf5Compression, outPath))

def runHalCons(halPath, outputPath):
    return runExecutable("halCons %s > outputPath" % halPath)

def initOptions():
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('--string', type=str,
                        help='string')
    parser.add_argument('--int', type=int, default=50, 
                        help='int. default=%(default)s')
    return parser
            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description='Get resources of executable')
    parser.add_argument('--string', type=str,
                        help='string')
    parser.add_argument('--int', type=int, default=50, 
                        help='int. default=%(default)s')
    parser = initOptions()
    args = parser.parse_args()
    rval = 0

    try:
        tempDir = getTempDirectory(rootDir="./")
        tempFile = getTempFile(rootDir=tempDir)
    except:
        traceback.print_exc(file=sys.stdout)
        exit(1)
    try:
        print runHalGen("large", 10000000, 5, tempFile)
        print runHalCons(tempFile, getTempFile(rootDir=tempDir))
    except:
        traceback.print_exc(file=sys.stdout)
        rval = 1

    system("rm -rf %s" % tempDir)
    return rval

if __name__ == "__main__":
    sys.exit(main())
