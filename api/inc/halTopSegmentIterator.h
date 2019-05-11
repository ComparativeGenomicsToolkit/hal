/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALTOPSEGMENTITERATOR_H
#define _HALTOPSEGMENTITERATOR_H

#include "halDefs.h"
#include "halGenome.h"
#include "halSegmentIterator.h"
#include "halTopSegment.h"
#include <iostream>

namespace hal {

    /**
     * Interface for top segment iterator exposes the top segment
     * interface and some new methods for jumping around the genome.
     * Always hidden in smart pointers in the public interface.
     */
    class TopSegmentIterator : public virtual SegmentIterator {
      public:
        /* constructor */
        TopSegmentIterator(TopSegment *topSegment, hal_offset_t startOffset = 0, hal_offset_t endOffset = 0,
                           bool reversed = false)
            : SegmentIterator(startOffset, endOffset, reversed), _topSegment(topSegment) {
        }

        /* destructor */
        virtual ~TopSegmentIterator() {
        }

        /** Return a new copy of the iterator */
        TopSegmentIteratorPtr clone() const;

        /** Copy an input iterator.  More efficient than the above methods
         * as no new iterator needs to be allocated
         * @param ts Iterator to copy */
        void copy(const TopSegmentIteratorPtr &topSegIt);

        /** Move the iterator to the child of a given bottom segment
         * @param bs Bottom segment whose child will be moved to
         * @param child Index of child in bottom segment's genome */
        void toChild(const BottomSegmentIteratorPtr &botSegIt, hal_size_t child);

        /** Move the iterator to the child of a given bottom segment
         * @param bs Bottom segment whose child will be moved to
         * @param childGenome genome of child in bottom segment */
        void toChildG(const BottomSegmentIteratorPtr &botSegIt, const Genome *childGenome);

        /** Given a bottom segment, move to the top segment that contains
         * its start position.  The genome remains unchanged.  The iterator
         * will be sliced accordingly (reversed state also taken into account)
         * @param bs Bottom segment to parse up from */
        void toParseUp(const BottomSegmentIteratorPtr &botSegIt);

        /** Return a pointer to the current TopSegment.  NOTE: changes when iterator is modified.  */
        TopSegment *getTopSegment() {
            return _topSegment.get();
        }

        /** Return a pointer to the current TopSegment.  NOTE: changes when iterator is modified.  */
        const TopSegment *getTopSegment() const {
            return _topSegment.get();
        }

        /** Return a pointer to the current TopSegment (terse).   NOTE: changes when iterator is modified.  */
        TopSegment *tseg() {
            return _topSegment.get();
        }

        /** Return a pointer to the current TopSegment (terse).   NOTE: changes when iterator is modified.  */
        const TopSegment *tseg() const {
            return _topSegment.get();
        }

        /** Test equality with other iterator (current implementation does not
         * take into account reverse state or offsets -- too review)
         * FIXME merge with operator==??
         * @param other Iterator to test equality to */
        bool equals(const TopSegmentIteratorPtr &other) const {
            assert(getGenome() == other->getGenome());
            return getArrayIndex() == other->getArrayIndex();
        }

        /* equality operator */
        bool operator==(const TopSegmentIterator &other) const {
            assert(_topSegment->getGenome() == other.getTopSegment()->getGenome());
            return getArrayIndex() == other.getArrayIndex();
        }

        /* inequality operator */
        bool operator!=(const TopSegmentIterator &other) const {
            return !(*this == other);
        }

        /** Move iterator to next paralgous segment.  Iterator will be reversed
        * if the next segment is in a different orientation wrt their common
        * parent */
        void toNextParalogy();

        // FIXME: document or change way getting segment works
        virtual Segment *getSegment() {
            return _topSegment.get();
        }
        virtual const Segment *getSegment() const {
            return _topSegment.get();
        }

        virtual void print(std::ostream &os) const;

      private:
        hal_size_t getNumSegmentsInGenome() const {
            return getGenome()->getNumTopSegments();
        }

        TopSegmentPtr _topSegment;
    };

    inline bool operator==(TopSegmentIteratorPtr p1, TopSegmentIteratorPtr p2) {
        if (p1.get() == NULL || p2.get() == NULL) {
            return p1.get() == NULL && p2.get() == NULL;
        }
        return p1->equals(p2);
    }

    inline bool operator!=(TopSegmentIteratorPtr p1, TopSegmentIteratorPtr p2) {
        return !(p1 == p2);
    }
}

#endif
// Local Variables:
// mode: c++
// End:
