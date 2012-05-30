#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt


import argparse
import os
import sys
import traceback
import time
import resource
import psutil

from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import popenCatch
from sonLib.bioio import system
    
def runHalGen(preset, hdf5Chunk, hdf5Compression, outPath):
    system("halRandGen --preset %s --hdf5Chunk %d\
    --hdf5Compression %d %s" % (preset, hdf5Chunk, hdf5Compression, outPath))

def runHalCons(halPath, outputPath):
    system("halCons %s > outputPath" % halPath)

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
        for chunkSize in [100000, 1000000, 10000000]:
            for compression in [0, 5, 9]:
                try:
                    tempDir = getTempDirectory(rootDir="./")
                    tempFile = getTempFile(suffix=".h5", rootDir=tempDir)
                except:
                    traceback.print_exc(file=sys.stdout)
                    exit(1)

                t = time.time()
                runHalGen("big", chunkSize, 5, tempFile)
                th = time.time() - t
                runHalCons(tempFile, getTempFile(rootDir=tempDir))
                tc = time.time() - th - t
                print "chunk=%d comp=%d:  generate %f.3  cons %f.3" % (
                    chunkSize, compression, th, tc)


    except:
        traceback.print_exc(file=sys.stdout)
        rval = 1

    system("rm -rf %s" % tempDir)
    return rval

if __name__ == "__main__":
    sys.exit(main())
