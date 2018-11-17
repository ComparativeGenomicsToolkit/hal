/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#ifndef _HALMAPPEDSEGMENTCONTAINERS_H
#define _HALMAPPEDSEGMENTCONTAINERS_H
#include "halDefs.h"
#include <set>

namespace hal {
    /* Functor for set compare; implemented in halMappedSegment.cpp This needs
     * to be in it's is used by hal::Segment which is required by
     * halMappedSegment.h.  It also needs to have the implementation split
     * from the class definition since the compare function references
     * MappedSegment. */
    struct MappedSegmentLess {
        bool operator()(const hal::MappedSegment& m1,
                        const hal::MappedSegment& m2) const;
        bool operator()(const hal::MappedSegmentPtr& m1,
                        const hal::MappedSegmentPtr& m2) const;
    };

    /* set of MappedSegments objects */
    class MappedSegmentSet: public std::set<MappedSegmentPtr, MappedSegmentLess> {
    };
}  

#endif

// Local Variables:
// mode: c++
// End:
