/*
 * Copyright (C) 2013 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HALMAPPEDSEGMENT_H
#define _HALMAPPEDSEGMENT_H

#include "halBottomSegmentIterator.h"
#include "halDefs.h"
#include "halSegmentIterator.h"
#include "halSlicedSegment.h"
#include "halTopSegmentIterator.h"
#include <list>

namespace hal {
    /**
     * Interface for a mapped segment.  A mapped segment stores a source segment
     * and a homologous region in a target genome (to which it was mapped).  Mapped
     * segments are used to keep pairwise alignment fragments across the tree as
     * an alternative to the column iterator.
     */
    class MappedSegment {
      public:
        /* Constructor */
        MappedSegment(SegmentIteratorPtr sourceSegIt, SegmentIteratorPtr targetSegIt);

        /** Destructor */
        virtual ~MappedSegment() {
        }

        /** Get the original segment from which this segment was mapped */
        SlicedSegment *getSource() const {
            return _source.get();
        }

        /** Get the original segment from which this segment was mapped */
        SlicedSegmentPtr getSourcePtr() const {
            return _source;
        }

        /** Get the original segment from which this segment was mapped */
        SegmentIteratorPtr getSourceIteratorPtr() {
            return _source;
        }

        /** get the target object to which this segment was mapped */
        SlicedSegment *getTarget() {
            return _target.get();
        }

        /** get the target object to which this segment was mapped */
        SlicedSegmentPtr getTargetPtr() {
            return _target;
        }

        /** get the target object to which this segment was mapped */
        SegmentIterator *getTargetIterator() {
            return _target.get();
        }

        /** get the target object to which this segment was mapped */
        SegmentIteratorPtr getTargetIteratorPtr() {
            return _target;
        }

        /** Comparison used to store in stl sets and maps.  We sort based
         * on the coordinate of the mapped segemnt's (target) interval as the primary
         * index and the target genome as the secondary index.  */
        virtual bool lessThan(const MappedSegment *other) const;

        /** Comparison used to store in STL containers.  We sort based
         * on the coordinate of the *Source* genome interval as the primary
         * index and the target genome as the secondary index.  */
        virtual bool lessThanBySource(const MappedSegment *other) const;

        /** Comparison used to determine uniqueness in lists.  Tests lessThan
         * in both directions.  Note that equality is the same regardless of
         * whether or not we use the source segment as our primary index. */
        virtual bool equals(const MappedSegment *other) const;

        /** Flip the mapping direction.  This segment becomes the source, and
         * the source becomes this.*/
        virtual void flip();

        /** Reverse both segments.  Also swap their start and end offsets.
         * Note that toReverse() does not reverse the source segment.*/
        virtual void fullReverse();

        /** Return of a copy of the mapped segment */
        virtual MappedSegment *clone() const;

        /** Test if mapped segment can be merged to the right with input
         * segment.  will return false if the right coordinate of this is in
         * either (optional) cutSet.*/
        virtual bool canMergeRightWith(const MappedSegmentPtr &nextSeg, const std::set<hal_index_t> *cutSet = NULL,
                                       const std::set<hal_index_t> *sourceCutSet = NULL) const;

        // FIXME: are non-ptr functors needed?
        /** Functor for sorted STL containers, sorting by origin as primary
         * index */
        struct LessSource {
            bool operator()(const MappedSegment &ms1, const MappedSegment &ms2) const {
                return ms1.lessThanBySource(&ms2);
            }
        };

        /** Functor for sorted STL containers, sorting by origin as primary
         * index */
        struct LessSourcePtr {
            bool operator()(const MappedSegmentPtr &ms1, const MappedSegmentPtr &ms2) const {
                return ms1->lessThanBySource(ms2.get());
            }
        };

        /** Functor for sorted STL containers, sorting  by target as primary
         * index */
        struct Less {
            bool operator()(const MappedSegment &ms1, const MappedSegment &ms2) const {
                return ms1.lessThan(&ms2);
            }
        };

        /** Functor for sorted STL containers, sorting  by target as primary
         * index */
        struct LessPtr {
            bool operator()(const MappedSegmentPtr &ms1, const MappedSegmentPtr &ms2) const {
                return ms1->lessThan(ms2.get());
            }
        };

        /** Functor for STL sorted lists to test for uniqueness */
        struct EqualTo {
            bool operator()(const MappedSegment &ms1, const MappedSegment &ms2) const {
                return ms1.equals(&ms2);
            }
        };

        /** Functor for STL sorted lists to test for uniqueness */
        struct EqualToPtr {
            bool operator()(const MappedSegmentPtr &ms1, const MappedSegmentPtr &ms2) const {
                return ms1->equals(ms2.get());
            }
        };

        // NEEDS TO BE ADDED TO SEGMENT INTERFACE
        virtual std::ostream &print(std::ostream &os) const;

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
        virtual bool isMissingData(double nThreshold) const;
        virtual bool isTop() const;

        // SLICED SEGMENT INTERFACE
        virtual void toReverse();
        virtual void toReverseInPlace();
        virtual hal_offset_t getStartOffset() const;
        virtual hal_offset_t getEndOffset() const;
        virtual void slice(hal_offset_t startOffset, hal_offset_t endOffset);
        virtual bool getReversed() const;

        // used by halSegmentMapper (FIXME: lef tover from transition, might want to clean up or doc)
        TopSegmentIteratorPtr targetAsTop() const;
        BottomSegmentIteratorPtr targetAsBottom() const;
        TopSegmentIteratorPtr sourceAsTop() const;
        BottomSegmentIteratorPtr sourceAsBottom() const;
        SegmentIteratorPtr sourceClone() const;
        void setTarget(SegmentIteratorPtr target) {
            _target = target;
        }

      private:
        static int fastComp(const SegmentIteratorPtr &s1, const SegmentIteratorPtr &s2);

        static int boundComp(const SegmentIteratorPtr &s1, const SegmentIteratorPtr &s2);

        static int slowComp(const SegmentIteratorPtr &s1, const SegmentIteratorPtr &s2);

      private:
        SegmentIteratorPtr _source;
        SegmentIteratorPtr _target;
    };
}

inline hal::TopSegmentIteratorPtr hal::MappedSegment::targetAsTop() const {
    return std::dynamic_pointer_cast<TopSegmentIterator>(_target);
}

inline hal::BottomSegmentIteratorPtr hal::MappedSegment::targetAsBottom() const {
    return std::dynamic_pointer_cast<BottomSegmentIterator>(_target);
}

inline hal::TopSegmentIteratorPtr hal::MappedSegment::sourceAsTop() const {
    return std::dynamic_pointer_cast<TopSegmentIterator>(_source);
}

inline hal::BottomSegmentIteratorPtr hal::MappedSegment::sourceAsBottom() const {
    return std::dynamic_pointer_cast<BottomSegmentIterator>(_source);
}

inline hal::SegmentIteratorPtr hal::MappedSegment::sourceClone() const {
    if (_source->isTop()) {
        return sourceAsTop()->clone();
    } else {
        return sourceAsBottom()->clone();
    }
}

#endif
// Local Variables:
// mode: c++
// End:
