#!/usr/bin/env python3

#Copyright (C) 2012 by Glenn Hickey
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt


import os
import sys
import time

from sonLib.bioio import popenCatch
from sonLib.bioio import getTotalCpuTimeAndMemoryUsage
from sonLib.bioio import system
            
def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) != 2:
        print("usage: runAndGetResources.py \'cmdline\'")
        exit(1)
    cmdline = argv[1]
    wallStart = time.time()
    output = popenCatch(cmdline)
    wallClock = time.time() - wallStart
    print((wallClock,) + getTotalCpuTimeAndMemoryUsage())
    return 0

if __name__ == "__main__":
    sys.exit(main())
