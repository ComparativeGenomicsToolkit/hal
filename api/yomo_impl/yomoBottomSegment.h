/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _YOMOBOTTOMSEGMENT_H
#define _YOMOBOTTOMSEGMENT_H

#include "halBottomSegment.h"
#include "halDefs.h"
#include "yomoExternalArray.h"
#include "yomoGenome.h"
#include <H5Cpp.h>

namespace hal {

    class YomoBottomSegment : public BottomSegment {
      public:
        /** Constructor
        * @param genome Smart pointer to genome to which segment belongs
        * @param array YOMO array containg segment
        * @param index Index of segment in the array */
        YomoBottomSegment(YomoGenome *genome, YomoExternalArray *array, hal_index_t index);

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

        // BOTTOM SEGMENT INTERFACE
        hal_size_t getNumChildren() const;
        hal_index_t getChildIndex(hal_size_t i) const;
        hal_index_t getChildIndexG(const Genome *childGenome) const;
        bool hasChild(hal_size_t child) const;
        bool hasChildG(const Genome *childGenome) const;
        void setChildIndex(hal_size_t i, hal_index_t childIndex);
        bool getChildReversed(hal_size_t i) const;
        void setChildReversed(hal_size_t child, bool isReversed);
        hal_index_t getTopParseIndex() const;
        void setTopParseIndex(hal_index_t parseIndex);
        hal_offset_t getTopParseOffset() const;
        bool hasParseUp() const;
        hal_index_t getLeftChildIndex(hal_size_t i) const;
        hal_index_t getRightChildIndex(hal_size_t i) const;

        // YOMO SPECIFIC
        static H5::CompType dataType(hal_size_t numChildren);
        static hal_size_t numChildrenFromDataType(const H5::DataType &dataType);

      private:
        YomoGenome *getYomoGenome() const {
            return static_cast<YomoGenome *>(_genome);
        }

        static const size_t genomeIndexOffset;
        static const size_t lengthOffset;
        static const size_t topIndexOffset;
        static const size_t firstChildOffset;
        static const size_t totalSize(hal_size_t numChildren);

        YomoExternalArray *_array;
    };

    // INLINE members
    inline void YomoBottomSegment::setArrayIndex(Genome *genome, hal_index_t arrayIndex) {
        _genome = dynamic_cast<YomoGenome *>(genome);
        assert(_genome != NULL);
        _array = &dynamic_cast<YomoGenome *>(_genome)->_bottomArray;
        assert(arrayIndex < (hal_index_t)_array->getSize());
        _index = arrayIndex;
    }

    inline hal_index_t YomoBottomSegment::getStartPosition() const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, genomeIndexOffset);
    }

    inline hal_index_t YomoBottomSegment::getEndPosition() const {
        assert(_index >= 0);
        return getStartPosition() + (hal_index_t)(getLength() - 1);
    }

    inline hal_size_t YomoBottomSegment::getLength() const {
        assert(_index >= 0);
        return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) -
               _array->getValue<hal_size_t>(_index, genomeIndexOffset);
    }

    inline const Sequence *YomoBottomSegment::getSequence() const {
        return _genome->getSequenceBySite(getStartPosition());
    }

    inline hal_size_t YomoBottomSegment::getNumChildren() const {
        return _genome->getNumChildren();
    }

    inline hal_index_t YomoBottomSegment::getChildIndex(hal_size_t i) const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)));
    }

    inline hal_index_t YomoBottomSegment::getChildIndexG(const Genome *childGenome) const {
        assert(_index >= 0);
        return getChildIndex(_genome->getChildIndex(childGenome));
    }

    inline bool YomoBottomSegment::hasChild(hal_size_t i) const {
        return getChildIndex(i) != NULL_INDEX;
    }

    inline bool YomoBottomSegment::hasChildG(const Genome *childGenome) const {
        return getChildIndexG(childGenome) != NULL_INDEX;
    }

    inline void YomoBottomSegment::setChildIndex(hal_size_t i, hal_index_t childIndex) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)), childIndex);
    }

    inline bool YomoBottomSegment::getChildReversed(hal_size_t i) const {
        assert(_index >= 0);
        return _array->getValue<bool>((hsize_t)_index,
                                      firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)) + sizeof(hal_index_t));
    }

    inline void YomoBottomSegment::setChildReversed(hal_size_t i, bool isReversed) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)) + sizeof(hal_index_t),
                         isReversed);
    }

    inline hal_index_t YomoBottomSegment::getTopParseIndex() const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, topIndexOffset);
    }

    inline void YomoBottomSegment::setTopParseIndex(hal_index_t parseIndex) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, topIndexOffset, parseIndex);
    }

    inline bool YomoBottomSegment::hasParseUp() const {
        return getTopParseIndex() != NULL_INDEX;
    }

    inline hal_index_t YomoBottomSegment::getArrayIndex() const {
        return _index;
    }

    inline bool YomoBottomSegment::isFirst() const {
        assert(getSequence() != NULL);
        return _index == 0 || _index == (hal_index_t)getSequence()->getBottomSegmentArrayIndex();
    }

    inline bool YomoBottomSegment::isLast() const {
        assert(getSequence() != NULL);
        return _index == (hal_index_t)_array->getSize() - 1 ||
               _index == getSequence()->getBottomSegmentArrayIndex() + (hal_index_t)getSequence()->getNumBottomSegments() - 1;
    }

    inline hal_index_t YomoBottomSegment::getLeftChildIndex(hal_size_t i) const {
        assert(isFirst() == false);
        YomoBottomSegment leftSeg(getYomoGenome(), _array, _index - 1);
        return leftSeg.getChildIndex(i);
    }

    inline hal_index_t YomoBottomSegment::getRightChildIndex(hal_size_t i) const {
        assert(isLast() == false);
        YomoBottomSegment rightSeg(dynamic_cast<YomoGenome *>(_genome), _array, _index + 1);
        return rightSeg.getChildIndex(i);
    }
}

#endif
// Local Variables:
// mode: c++
// End:
