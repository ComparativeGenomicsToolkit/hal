/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALGAPPEDTOPSEGMENTITERATOR_H
#define _HALGAPPEDTOPSEGMENTITERATOR_H

#include "halDefs.h"
#include "halGappedSegmentIterator.h"
#include "halTopSegmentIterator.h"
#include <iostream>

namespace hal {

    /**
     * Interface for Gapped Top Segment iterator.  Only used internally
     * and probably shouldn't be in public interface.
     */
    class GappedTopSegmentIterator : virtual public GappedSegmentIterator {
      public:
        /** constructor */
        GappedTopSegmentIterator(TopSegmentIteratorPtr leftTopSegIt, hal_size_t gapThreshold, bool atomic);

        /** Destructor */
        virtual ~GappedTopSegmentIterator() {
        }

        /** Return a copy of the iterator */
        virtual GappedTopSegmentIteratorPtr clone() const;

        /** Copy another iterator into the current iterator (more efficient
         * than above methods since no new iterators are created */
        virtual void copy(const GappedTopSegmentIteratorPtr &gapTopSegIt);

        /** Move to child of given iterator
         * @param bs Move to child of this iterator */
        virtual void toChild(const GappedBottomSegmentIteratorPtr &gapBotSegIt);

        /** Test equality with other iterator
         * @param other */
        virtual bool equals(const GappedTopSegmentIteratorPtr &other) const;

        /** Test if iterator abuts other iterator */
        virtual bool adjacentTo(const GappedTopSegmentIteratorPtr &other) const;

        /** Test if iterator has a parent */
        virtual bool hasParent() const;

        /** Test if iterator has a duplicate */
        virtual bool hasNextParalogy() const;

        /** Move to next paralogy */
        virtual void toNextParalogy() const;

        /** Test if in reverse orientationt with respect to parent */
        virtual bool getParentReversed() const;

        /** Return the leftmost segment of the iterator
         * (note that moving the returned iterator will corrupt the
         * current gapped iterator.  this is a bug) */
        virtual TopSegmentIteratorPtr getLeft() const;

        /** Return the rightmost segment of the iterator
         * (note that moving the returned iterator will corrupt the
         * current gapped iterator.  this is a bug) */
        virtual TopSegmentIteratorPtr getRight() const;

        /** Reset the gapped iterator.
         * @param ts This will be the left segment of the current iterator. The
         * right segment will be extended as far as possible */
        virtual void setLeft(const TopSegmentIteratorPtr &topSegIt);

        /** For every set of paralogous top segments in a given genome, we identify a
         * unique segment as the canonical reference.  This is the segment that
         * will be traversed when disabling duplication edges when mapping via
         * the column iterator or mappedSegment interface.  It is also the segment
         * that is connected from its parent's down edge.*/
        virtual bool isCanonicalParalog() const;

        /* get the current segment */
        virtual Segment *getSegment();

        /* get the current segment */
        virtual const Segment *getSegment() const;

        // SEGMENT INTERFACE
        virtual void setArrayIndex(Genome *genome, hal_index_t arrayIndex);
        virtual const Genome *getGenome() const;
        virtual Genome *getGenome();
        virtual const Sequence *getSequence() const;
        virtual hal_index_t getStartPosition() const;
        virtual hal_index_t getEndPosition() const;
        virtual hal_size_t getLength() const;
        virtual void getString(std::string &outString) const;
        virtual void setCoordinates(hal_index_t startPos, hal_size_t length);
        virtual hal_index_t getArrayIndex() const;
        virtual bool leftOf(hal_index_t genomePos) const;
        virtual bool rightOf(hal_index_t genomePos) const;
        virtual bool overlaps(hal_index_t genomePos) const;
        virtual bool isFirst() const;
        virtual bool isLast() const;
        virtual bool isTop() const;
        virtual std::ostream &print(std::ostream &os) const;

        // SEGMENT ITERATOR INTERFACE
        virtual void toLeft(hal_index_t leftCutoff = NULL_INDEX);
        virtual void toRight(hal_index_t rightCutoff = NULL_INDEX);
        virtual void toReverse();
        virtual void toReverseInPlace();
        virtual void toSite(hal_index_t position, bool slice = true);
        virtual hal_offset_t getStartOffset() const;
        virtual hal_offset_t getEndOffset() const;
        virtual void slice(hal_offset_t startOffset, hal_offset_t endOffset);
        virtual bool getReversed() const;

        // GAPPED SEGMENT ITERATOR INTERFACE
        virtual hal_size_t getGapThreshold() const;
        virtual bool getAtomic() const;
        virtual hal_size_t getChildIndex() const;
        virtual hal_size_t getNumSegments() const;
        virtual hal_size_t getNumGaps() const;
        virtual hal_size_t getNumGapBases() const;
        virtual hal_index_t getLeftArrayIndex() const;
        virtual hal_index_t getRightArrayIndex() const;

      private:
        bool compatible(TopSegmentIteratorPtr leftTopSegIt, TopSegmentIteratorPtr rightTopSegIt) const;

        void extendRight();
        void extendLeft();

        void toLeftNextUngapped(BottomSegmentIteratorPtr botSegIt) const;
        void toRightNextUngapped(BottomSegmentIteratorPtr botSegIt) const;
        void toLeftNextUngapped(TopSegmentIteratorPtr topSeqIt) const;
        void toRightNextUngapped(TopSegmentIteratorPtr topSeqIt) const;

        // keep convention of other iterators where const-ness only applies
        // to the database and not the iterator...
        TopSegmentIteratorPtr _left;
        TopSegmentIteratorPtr _right;
        BottomSegmentIteratorPtr _leftParent;
        BottomSegmentIteratorPtr _rightParent;
        TopSegmentIteratorPtr _leftDup;
        TopSegmentIteratorPtr _rightDup;
        TopSegmentIteratorPtr _temp;
        TopSegmentIteratorPtr _temp2;
        hal_size_t _childIndex;
        hal_size_t _gapThreshold;
        bool _atomic;
    };

    inline bool operator==(GappedTopSegmentIteratorPtr p1, GappedTopSegmentIteratorPtr p2) {
        if (p1.get() == NULL || p2.get() == NULL) {
            return p1.get() == NULL && p2.get() == NULL;
        }
        return p1->equals(p2);
    }

    inline bool operator!=(GappedTopSegmentIteratorPtr p1, GappedTopSegmentIteratorPtr p2) {
        return !(p1 == p2);
    }
}

#endif
// Local Variables:
// mode: c++
// End:
