/*
 * HAL format 3: conversion from a HAL (HDF5 or mmap) alignment.
 *
 * This is the only part of format 3 that uses the HAL API, and so the only part that
 * links HDF5.  Everything else (hal3Manifest, hal3Validate) depends on liftmap and
 * twobit64 alone.
 */
#ifndef _HAL3CONVERT_H
#define _HAL3CONVERT_H

#include "halAlignment.h"
#include <cstdint>
#include <ostream>
#include <string>

namespace hal3 {

struct ConvertOptions {
    uint32_t chunkRuns = 8192;          /* runs per chunk at most */
    uint64_t chunkParentSpan = 1000000; /* parent-axis span per chunk at most */
    bool dna = true;                    /* write the .2bit (when the genome has DNA) */
};

struct GenomeConversion {
    uint64_t topSegments = 0;    /* in HAL */
    uint64_t alignedSegments = 0; /* of those, with a parent */
    uint64_t alignedBases = 0;
    uint64_t runs = 0;           /* maximal runs written */
    bool wroteDna = false;
};

/* Write <genomesDir>/<stem>.hgn, and <stem>.2bit if the genome has DNA and
 * opts.dna.  Throws hal_exception or std::runtime_error. */
GenomeConversion convertGenome(const hal::Alignment *alignment, const std::string &genomeName,
                               const std::string &genomesDir, const ConvertOptions &opts);

/* Compare the written files with HAL: the sequence tables, then nSamples random
 * bases, each checked for its parent position and strand (or for having none) and
 * for its DNA.  Returns the number of mismatches, and describes the first few on
 * log. */
uint64_t checkGenome(const hal::Alignment *alignment, const std::string &genomeName,
                     const std::string &genomesDir, uint64_t nSamples, uint64_t seed, std::ostream &log);
}

#endif
