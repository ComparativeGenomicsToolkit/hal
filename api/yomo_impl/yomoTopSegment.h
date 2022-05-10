/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOTOPSEGMENT_H
#define _YOMOTOPSEGMENT_H

#include "halTopSegment.h"
#include "yomoExternalArray.h"
#include "yomoGenome.h"
#include <H5Cpp.h>

namespace hal {

    class YomoTopSegment : public TopSegment {
      public:
        /** Constructor
         * @param genome Smart pointer to genome to which segment belongs
         * @param array YOMO array containg segment
         * @param index Index of segment in the array */
        YomoTopSegment(YomoGenome *genome, YomoExternalArray *array, hal_index_t index);

        // SEGMENT INTERFACE
        void setArrayIndex(Genome *genome, hal_index_t arrayIndex);
        const Sequence *getSequence() const;
        hal_index_t getStartPosition() const;
        hal_index_t getEndPosition() const;
        hal_size_t getLength() const;
        void setCoordinates(hal_index_t startPos, hal_size_t length);
        hal_index_t getArrayIndex() const;
        bool isFirst() const;
        bool isLast() const;
        void print(std::ostream &os) const;

        // TOP SEGMENT INTERFACE
        hal_index_t getParentIndex() const;
        bool hasParent() const;
        void setParentIndex(hal_index_t parIdx);
        bool getParentReversed() const;
        void setParentReversed(bool isReversed);
        hal_index_t getBottomParseIndex() const;
        void setBottomParseIndex(hal_index_t botParseIdx);
        hal_offset_t getBottomParseOffset() const;
        bool hasParseDown() const;
        hal_index_t getNextParalogyIndex() const;
        bool hasNextParalogy() const;
        void setNextParalogyIndex(hal_index_t parIdx);
        hal_index_t getLeftParentIndex() const;
        hal_index_t getRightParentIndex() const;
        bool isCanonicalParalog() const;

        // YOMO SPECIFIC
        static H5::CompType dataType();

      private:
        YomoGenome *getYomoGenome() const {
            return static_cast<YomoGenome *>(_genome);
        }

        static const size_t genomeIndexOffset;
        static const size_t bottomIndexOffset;
        static const size_t parIndexOffset;
        static const size_t parentIndexOffset;
        static const size_t parentReversedOffset;
        static const size_t totalSize;

        YomoExternalArray *_array;
    };

    // INLINE members
    inline void YomoTopSegment::setArrayIndex(Genome *genome, hal_index_t arrayIndex) {
        _genome = dynamic_cast<YomoGenome *>(genome);
        assert(_genome != NULL);
        _array = &dynamic_cast<YomoGenome *>(_genome)->_topArray;
        assert(arrayIndex < (hal_index_t)_array->getSize());
        _index = arrayIndex;
    }

    inline hal_index_t YomoTopSegment::getStartPosition() const {
        return _array->getValue<hal_index_t>(_index, genomeIndexOffset);
    }

    inline hal_index_t YomoTopSegment::getEndPosition() const {
        return getStartPosition() + (hal_index_t)(getLength() - 1);
    }

    inline hal_size_t YomoTopSegment::getLength() const {
        return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) -
               _array->getValue<hal_size_t>(_index, genomeIndexOffset);
    }

    inline const Sequence *YomoTopSegment::getSequence() const {
        return _genome->getSequenceBySite(getStartPosition());
    }

    inline bool YomoTopSegment::hasParseDown() const {
        return getBottomParseIndex() != NULL_INDEX;
    }

    inline hal_index_t YomoTopSegment::getNextParalogyIndex() const {
        return _array->getValue<hal_index_t>(_index, parIndexOffset);
    }

    inline bool YomoTopSegment::hasNextParalogy() const {
        return getNextParalogyIndex() != NULL_INDEX;
    }

    inline void YomoTopSegment::setNextParalogyIndex(hal_index_t parIdx) {
        assert(parIdx != _index);
        _array->setValue(_index, parIndexOffset, parIdx);
    }

    inline hal_index_t YomoTopSegment::getParentIndex() const {
        return _array->getValue<hal_index_t>(_index, parentIndexOffset);
    }

    inline bool YomoTopSegment::hasParent() const {
        return getParentIndex() != NULL_INDEX;
    }

    inline void YomoTopSegment::setParentIndex(hal_index_t parentIndex) {
        _array->setValue(_index, parentIndexOffset, parentIndex);
    }

    inline bool YomoTopSegment::getParentReversed() const {
        return _array->getValue<bool>(_index, parentReversedOffset);
    }

    inline void YomoTopSegment::setParentReversed(bool isReversed) {
        _array->setValue(_index, parentReversedOffset, isReversed);
    }

    inline hal_index_t YomoTopSegment::getBottomParseIndex() const {
        return _array->getValue<hal_index_t>(_index, bottomIndexOffset);
    }

    inline void YomoTopSegment::setBottomParseIndex(hal_index_t parseIndex) {
        _array->setValue(_index, bottomIndexOffset, parseIndex);
    }

    inline hal_index_t YomoTopSegment::getArrayIndex() const {
        return _index;
    }

    inline bool YomoTopSegment::isFirst() const {
        assert(getSequence() != NULL);
        return _index == 0 || _index == (hal_index_t)getSequence()->getTopSegmentArrayIndex();
    }

    inline bool YomoTopSegment::isLast() const {
        assert(getSequence() != NULL);
        return _index == (hal_index_t)_array->getSize() - 1 ||
               _index == getSequence()->getTopSegmentArrayIndex() + (hal_index_t)getSequence()->getNumTopSegments() - 1;
    }

    inline hal_index_t YomoTopSegment::getLeftParentIndex() const {
        assert(isFirst() == false);
        YomoTopSegment leftSeg(dynamic_cast<YomoGenome *>(_genome), _array, _index - 1);
        return leftSeg.getParentIndex();
    }

    inline hal_index_t YomoTopSegment::getRightParentIndex() const {
        assert(isLast() == false);
        YomoTopSegment rightSeg(dynamic_cast<YomoGenome *>(_genome), _array, _index + 1);
        return rightSeg.getParentIndex();
    }
}

#endif
// Local Variables:
// mode: c++
// End:
