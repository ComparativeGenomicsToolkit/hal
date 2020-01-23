#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt


import argparse
import os
import sys
import traceback
import time
import random
import resource
import psutil

from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import popenCatch
from sonLib.bioio import system
    
def runHalGen(preset, seed, hdf5Chunk, hdf5Compression, outPath):
    system("halRandGen --preset %s --seed %d --hdf5Chunk %d\
    --hdf5Compression %d %s" % (preset, seed, hdf5Chunk, hdf5Compression, outPath))

def runHalCons(halPath, outputPath):
    system("halCons %s > outputPath" % halPath)

            
def main(argv=None):
    if argv is None:
        argv = sys.argv

    seed = random.randint(0, 2**31)
    parser = argparse.ArgumentParser(description='Run little hal test')
    parser.add_argument('--preset', type=str,
                        help='halGenRandom preset to use [small, medium, big, large]', default='small')
    args = parser.parse_args()
    rval = 0
    print("chunk, comp, time(gen), time(cons), fsize(k)")
    try:
        for chunkSize in [10000, 100000, 1000000, 10000000]:
            for compression in [0, 2, 5, 7, 9]:
                try:
                    tempDir = getTempDirectory(rootDir="./")
                    tempFile = getTempFile(suffix=".h5", rootDir=tempDir)
                except:
                    traceback.print_exc(file=sys.stdout)
                    return 1

                t = time.time()
                runHalGen(args.preset, seed, chunkSize, compression, tempFile)
                fsize = os.path.getsize(tempFile)
                th = time.time() - t
                runHalCons(tempFile, getTempFile(rootDir=tempDir))
                tc = time.time() - th - t
                print("%d, %d, %f.3, %f.3, %f.2" % (
                    chunkSize, compression, th, tc, fsize / 1024.))

    except:
        traceback.print_exc(file=sys.stdout)
        return 1

    system("rm -rf %s" % tempDir)
    return rval

if __name__ == "__main__":
    sys.exit(main())
