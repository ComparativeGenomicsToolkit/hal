#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""Take care of Level-of-detail files
"""

import os
from sonLib.bioio import system  
from optparse import OptionGroup

def fixLodFilePath(lodtxtfile, localHalfile, outdir):
    #fix the path of the original hal file to point to the created
    #link relative to the output directory
    relPath = os.path.relpath(localHalfile, start=outdir)
    lodTxtBuf = ''
    for line in open(lodtxtfile):
        tokens = line.split()
        if len(tokens) == 2 and tokens[0] == '0':
            lodTxtBuf += '0 %s\n' % relPath
        else:
            lodTxtBuf += line
    with open(lodtxtfile, 'w') as lodFile:
        lodFile.write(lodTxtBuf)
    
def getLodFiles(localHalfile, options, outdir):
    lodtxtfile = os.path.join(outdir, "lod.txt") #outdir/lod.txt
    loddir = os.path.join(outdir, "lod") #outdir/lod
    if options.lodtxtfile and options.loddir: #if lod files were given, then just make soft links to them
        if os.path.exists(lodtxtfile):
            if os.path.abspath(lodtxtfile) != os.path.abspath(options.lodtxtfile):
                system("rm %s" %lodtxtfile)
                system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.lodtxtfile), lodtxtfile))

        if os.path.exists(loddir):
            if os.path.abspath(loddir) != os.path.abspath(options.loddir):
                if os.path.islink(loddir):
                    system("rm %s" %loddir)
                else:
                    system("rm -Rf %s" %loddir)
                loddir = os.path.join(outdir, os.path.basename(options.loddir))
                system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
        else:
            system("ln -s %s %s" %(os.path.abspath(options.loddir), loddir))
    else: #if lod files were not given, create them using halLodInterpolate.py
        system("halLodInterpolate.py %s %s --outHalDir %s %s" %(localHalfile, lodtxtfile, loddir, options.lodOpts))
        fixLodFilePath(lodtxtfile, localHalfile, outdir)
    return lodtxtfile, loddir

def getLod(options, localHalfile, outdir):
    #Create lod files if useLod is specified
    lodtxtfile = ''
    loddir = ''
    options.lodOpts = ''
    if options.lodMaxBlock is not None:
        options.lodOpts += '--maxBlock %d ' % options.lodMaxBlock
    if options.lodScale is not None:
        options.lodOpts += '--scale %f ' % options.lodScale
    if options.lodMaxDNA is not None:
        options.lodOpts += '--maxDNA %d ' % options.lodMaxDNA
    if options.lodInMemory is True:
        options.lodOpts += '--inMemory '
    if options.lodMinSeqFrac is not None:
        options.lodOpts += '--minSeqFrac %f ' % options.lodMinSeqFrac
    if options.lodMinCovFrac is not None:
        options.lodOpts += '--minCovFrac %f ' % options.lodMinCovFrac
    if options.lodChunk is not None:
        options.lodOpts += '--chunk %d ' % options.lodChunk
    if options.maxCores and int(options.maxCores) > 1 and (options.lod or len(options.lodOpts) > 0):
        options.lodOpts += '--numProc %d ' % int(options.maxCores)
    if len(options.lodOpts) > 0:
        print(options.lodOpts)
        options.lod = True
    if options.lod:
        lodtxtfile, loddir = getLodFiles(localHalfile, options, outdir)
    return lodtxtfile, loddir

def getLodLowestLevel(lodtxtfile):
    """Gets the lowest level at which an LOD hal is used instead of the
    base-level hal file."""
    f = open(lodtxtfile, 'r')
    line = f.readline()
    level = int(line.split()[0])
    while line and level == 0:
        line = f.readline()
        if len(line.strip()) == 0:
            continue
        fields = line.strip().split()
        if len(fields) == 2 and fields[1] != "max":
            level = int(line.split()[0])
    f.close()
    return level

def addLodOptions(parser):
    group = parser.add_argument_group("LEVEL OF DETAILS", "Level-of-detail (LOD) options.")
    group.add_argument('--lod', dest='lod', action="store_true", default=False, help='If specified, create "level of detail" (lod) hal files and will put the lod.txt at the bigUrl instead of the original hal file. ')
    group.add_argument('--lodTxtFile', dest='lodtxtfile', help='"hal Level of detail" lod text file. If specified, will put this at the bigUrl instead of the hal file. ')
    group.add_argument('--lodDir', dest='loddir', help='"hal Level of detail" lod dir. If specified, will put this at the bigUrl instead of the hal file. ')
    group.add_argument('--lodMaxBlock', dest='lodMaxBlock', type=int, help='Maximum number of blocks to display in a hal level of detail (see halLodInterpolate.py --help for the default value).', default=None)
    group.add_argument('--lodScale', dest='lodScale', type=float, help='Scaling factor between two successive levels of detail (see halLodInterpolate.py --help for the default value).', default=None)
    group.add_argument('--lodMaxDNA', dest='lodMaxDNA', type=int, help='Maximum query length such that its hal level of detail will contain nucleotide information.  (see halLodInterpolate.py --help for the default value).', default=None)
    group.add_argument('--lodInMemory', dest='lodInMemory', action='store_true', help='Load entire hal file into memory when generating levels of detail instead of using hdf5 cache. Could result in drastic speedup. .', default=False)
    group.add_argument('--lodMinSeqFrac', dest='lodMinSeqFrac', type=float, help='Minumum sequence length to sample as fraction of step size for level of detail generation: ie sequences with length <= floor(minSeqFrac * step) are ignored (see halLodExtract --help for default value).', default=None)
    group.add_argument('--lodMinCovFrac', dest='lodMinCovFrac', type=float, help='Minimum fraction of a genome that must be covered by sequences that exceed --minSeqFrac * step.  LODs that would violate this threshold will not be generated (or displayed in  the browser).  This is seen a better than the alternative, which is to produce unreasonably sparse LODs because half the sequences were not sampled (see halLodInterpolate.py --help for default value).', default=None)
    group.add_argument('--lodChunk', dest='lodChunk', type=int, help='HDF5 chunk size for generated levels of detail (see halLodExtract --help for default value).', default=None)
    #group.add_argument('--snpwidth', dest='snpwidth', type=int, default=5000, help='Maximum window size to display SNPs. ')
    group = parser.add_argument_group(group)

