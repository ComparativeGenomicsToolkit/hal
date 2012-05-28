#!/usr/bin/env python

#Copyright (C) 2012 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt


import argparse
import os
import sys
import traceback
import time
import subprocess

from procWatch import monitor
from sonLib.bioio import getTempDirectory
from sonLib.bioio import getTempFile
from sonLib.bioio import system

def runExecutable(command):
    process = subprocess.Popen(command,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=sys.stderr,
                               bufsize=-1)

    class Arg:
        def __init__(self):
            pass
    procWatchArgs = Arg()
    procWatchArgs.pids = [process.pid]
    procWatchArgs.delay = 5.0
    procWatchArgs.cmds = [command]
    monitor(procWatchArgs)

    # below should be unecessary since monitor waits out process.
    # hopefully it can catch errors?
    sts = process.wait()

    if sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" %
                           (command, sts))
    
    return sts
    
def runHalGen(preset, hdf5Chunk, hdf5Compression, outPath):
    runExecutable("halRandGen --preset %s --hdf5Chunk %d\
    --hdf5Compression %d %s" % (preset, hdf5Chunk, hdf5Compression, outPath))

def runHalCons(halPath, outputPath):
    runExecutable("halCons %s > outputPath" % halPath)

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
    parser = initOptions()
    args = parser.parse_args()

    try:
        tempDir = getTempDirectory(rootDir="./")
        tempFile = getTempFile(rootDir=tempDir)
    except:
        traceback.print_exc(file=sys.stdout)
        exit(1)

    try:
        runHalGen("medium", 10000000, 5, tempFile)
        runHalCons(tempFile, getTempFile(rootDir=tempDir))
    except:
        traceback.print_exc(file=sys.stdout)
        exit(1)

    system("rm -rf %s" % tempDir)
    
if __name__ == "__main__":
    sys.exit(main())
