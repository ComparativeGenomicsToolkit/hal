/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALRANDOMDATA_H
#define _HALRANDOMDATA_H

#include <set>
#include "hal.h"

namespace hal {

void createRandomAlignment(Alignment* emptyAlignment,
                           double meanDegree,
                           double maxBranchLength,
                           hal_size_t maxGenomes,
                           hal_size_t minSegmentLength,
                           hal_size_t maxSegmentLength,
                           hal_size_t minSegments,
                           hal_size_t maxSegments,
                           int seed = -1);

void createRandomTree(Alignment* emptyAlignment,
                      double meanDegree,
                      double maxBranchLength,
                      hal_size_t maxGenomes);

void createRandomDimensions(Alignment* alignment,
                            hal_size_t minSegmentLength,
                            hal_size_t maxSegmentLength,
                            hal_size_t minSegments,
                            hal_size_t maxSegments);

void createRandomGenome(Alignment* alignment, Genome* genome);

void createRandomSegment(Genome* genome, 
                         hal_size_t indexInParent,
                         std::set<std::pair<hal_index_t, hal_index_t> >& 
                         edgeSet, 
                         TopSegmentIteratorPtr topIter, 
                         BottomSegmentIteratorPtr botIter,
                         double branchLength);

}
#endif
