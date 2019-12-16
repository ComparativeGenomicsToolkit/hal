#!/usr/bin/env python3

import os
import sys
import argparse

def getBedLength(bedPath):
    length = 0
    bedFile = open(bedPath)
    for line in bedFile:
        clnLine = line.strip()
        if len(clnLine) > 0 and clnLine[0] != "#":
            toks = clnLine.split()
            if len(toks) > 2:
                start = int(toks[1])
                end = int(toks[2])
                lineLength = end - start
                length += lineLength
    bedFile.close()
    return length

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser()
    parser.add_argument("bedPath", type=str,
                        help="path to bed file")
    
    args = parser.parse_args()
    print(getBedLength(args.bedPath))
        
if __name__ == "__main__":
    sys.exit(main())
