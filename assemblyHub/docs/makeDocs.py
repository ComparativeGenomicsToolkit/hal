#!/usr/bin/env python3

#Copyright (C) 2013 by Ngan Nguyen
# Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
#
#Released under the MIT license, see LICENSE.txt

"""
Make track documentation html files
"""
from hal.assemblyHub.docs.gcPercentDocs import makeGcPercentDocs 
from hal.assemblyHub.docs.alignabilityDocs import makeAlignabilityDocs
from hal.assemblyHub.docs.conservationDocs import makeConservationDocs
from hal.assemblyHub.docs.repeatMaskerDocs import makeRepeatMaskerDocs
from hal.assemblyHub.docs.hubCentralDocs import makeHubCentralDocs

def writeDocFiles(outdir, options):
    if options.gcContent:
        makeGcPercentDocs(outdir)
    if options.alignability:
        makeAlignabilityDocs(outdir)
    if options.conservation:
        makeConservationDocs(outdir)
    if options.rmskdir:
        makeRepeatMaskerDocs(outdir)
    makeHubCentralDocs(outdir)
    return

