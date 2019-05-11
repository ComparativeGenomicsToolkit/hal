/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALRANDOMDATA_H
#define _HALRANDOMDATA_H

#include "hal.h"
#include <set>

namespace hal {
    class RandNumberGen;

    void createRandomAlignment(RandNumberGen &rng, Alignment *emptyAlignment, double meanDegree, double maxBranchLength,
                               hal_size_t minGenomes, hal_size_t maxGenomes, hal_size_t minSegmentLength,
                               hal_size_t maxSegmentLength, hal_size_t minSegments, hal_size_t maxSegments);

    void createRandomTree(RandNumberGen &rng, Alignment *emptyAlignment, double meanDegree, double maxBranchLength,
                          hal_size_t minGenomes, hal_size_t maxGenomes);

    void createRandomDimensions(RandNumberGen &rng, Alignment *alignment, hal_size_t minSegmentLength,
                                hal_size_t maxSegmentLength, hal_size_t minSegments, hal_size_t maxSegments);

    void createRandomGenome(RandNumberGen &rng, Alignment *alignment, Genome *genome);

    void createRandomSegment(RandNumberGen &rng, Genome *genome, hal_size_t indexInParent,
                             std::set<std::pair<hal_index_t, hal_index_t>> &edgeSet, TopSegmentIteratorPtr topIter,
                             BottomSegmentIteratorPtr botIter, double branchLength);
}
#endif
// Local Variables:
// mode: c++
// End:
