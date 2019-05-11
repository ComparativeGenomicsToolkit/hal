/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDSEGMENTITERATOR_H
#define _HALGAPPEDSEGMENTITERATOR_H

#include "halSegmentIterator.h"

namespace hal {

    /**
     * Interface for general gappedSegment iterator.  Behaves like
     * a regular iterator, but operates on a linear sequence of
     * segments that are consistent modulo gaps.  Is only really used
     * internally at the moment so I haven't spent the effort to
     * make gapped iterators fully general (such as including the
     * Top/Bottom Segment iterfaces)
     */
    class GappedSegmentIterator : virtual public SegmentIterator {
      public:
        /** Destructor */
        virtual ~GappedSegmentIterator() {
        }

        /** Get the gap length threshold.  This is the maximum length (in sites)
        * of an indel such that it can be considered a gap (and therefore ignored
        * in the rearrangement analysis */
        virtual hal_size_t getGapThreshold() const = 0;

        /** When the gap iterator is in atomic mode, it allows for now gaps.
         * Also, it will not merge consistent segments together if they are adjacent
         * and there are no gaps (something that could happen when the threshold is
         * set to zero but atomic is false).  This is mostly a hack to get the gapped
         * iterator to behive as a single segment */
        virtual bool getAtomic() const = 0;

        /** Gapped iterators are tied to a specific parent child pair.  The child
         * index is the index of the child genome within the parent */
        virtual hal_size_t getChildIndex() const = 0;

        /** Get the number of segments that have been agglomerated together within
         * the gapped iterator */
        virtual hal_size_t getNumSegments() const = 0;

        /** Get the number of gaps within segments that have been agglomerated
         * together within the gapped iterator */
        virtual hal_size_t getNumGaps() const = 0;

        /** Get the number of bases within gaps within segments that have been
         * agglomerated together within the gapped iterator */
        virtual hal_size_t getNumGapBases() const = 0;

        /** Get the Segment array index of the left segment of the iterator */
        virtual hal_index_t getLeftArrayIndex() const = 0;

        /** Get the Segment array index of the right segment of the iterator */
        virtual hal_index_t getRightArrayIndex() const = 0;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
