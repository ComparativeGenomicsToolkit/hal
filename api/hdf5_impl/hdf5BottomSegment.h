/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5BOTTOMSEGMENT_H
#define _HDF5BOTTOMSEGMENT_H

#include "halBottomSegment.h"
#include "halDefs.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include <H5Cpp.h>

namespace hal {

    class Hdf5BottomSegment : public BottomSegment {
      public:
        /** Constructor
        * @param genome Smart pointer to genome to which segment belongs
        * @param array HDF5 array containg segment
        * @param index Index of segment in the array */
        Hdf5BottomSegment(Hdf5Genome *genome, Hdf5ExternalArray *array, hal_index_t index);

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

        // HDF5 SPECIFIC
        static H5::CompType dataType(hal_size_t numChildren);
        static hal_size_t numChildrenFromDataType(const H5::DataType &dataType);

      private:
        Hdf5Genome *getHdf5Genome() const {
            return static_cast<Hdf5Genome *>(_genome);
        }

        static const size_t genomeIndexOffset;
        static const size_t lengthOffset;
        static const size_t topIndexOffset;
        static const size_t firstChildOffset;
        static const size_t totalSize(hal_size_t numChildren);

        Hdf5ExternalArray *_array;
    };

    // INLINE members
    inline void Hdf5BottomSegment::setArrayIndex(Genome *genome, hal_index_t arrayIndex) {
        _genome = dynamic_cast<Hdf5Genome *>(genome);
        assert(_genome != NULL);
        _array = &dynamic_cast<Hdf5Genome *>(_genome)->_bottomArray;
        assert(arrayIndex < (hal_index_t)_array->getSize());
        _index = arrayIndex;
    }

    inline hal_index_t Hdf5BottomSegment::getStartPosition() const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, genomeIndexOffset);
    }

    inline hal_index_t Hdf5BottomSegment::getEndPosition() const {
        assert(_index >= 0);
        return getStartPosition() + (hal_index_t)(getLength() - 1);
    }

    inline hal_size_t Hdf5BottomSegment::getLength() const {
        assert(_index >= 0);
        return _array->getValue<hal_size_t>(_index + 1, genomeIndexOffset) -
               _array->getValue<hal_size_t>(_index, genomeIndexOffset);
    }

    inline const Sequence *Hdf5BottomSegment::getSequence() const {
        return _genome->getSequenceBySite(getStartPosition());
    }

    inline hal_size_t Hdf5BottomSegment::getNumChildren() const {
        return _genome->getNumChildren();
    }

    inline hal_index_t Hdf5BottomSegment::getChildIndex(hal_size_t i) const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)));
    }

    inline hal_index_t Hdf5BottomSegment::getChildIndexG(const Genome *childGenome) const {
        assert(_index >= 0);
        return getChildIndex(_genome->getChildIndex(childGenome));
    }

    inline bool Hdf5BottomSegment::hasChild(hal_size_t i) const {
        return getChildIndex(i) != NULL_INDEX;
    }

    inline bool Hdf5BottomSegment::hasChildG(const Genome *childGenome) const {
        return getChildIndexG(childGenome) != NULL_INDEX;
    }

    inline void Hdf5BottomSegment::setChildIndex(hal_size_t i, hal_index_t childIndex) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)), childIndex);
    }

    inline bool Hdf5BottomSegment::getChildReversed(hal_size_t i) const {
        assert(_index >= 0);
        return _array->getValue<bool>((hsize_t)_index,
                                      firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)) + sizeof(hal_index_t));
    }

    inline void Hdf5BottomSegment::setChildReversed(hal_size_t i, bool isReversed) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, firstChildOffset + i * (sizeof(hal_index_t) + sizeof(bool)) + sizeof(hal_index_t),
                         isReversed);
    }

    inline hal_index_t Hdf5BottomSegment::getTopParseIndex() const {
        assert(_index >= 0);
        return _array->getValue<hal_index_t>((hsize_t)_index, topIndexOffset);
    }

    inline void Hdf5BottomSegment::setTopParseIndex(hal_index_t parseIndex) {
        assert(_index >= 0);
        _array->setValue((hsize_t)_index, topIndexOffset, parseIndex);
    }

    inline bool Hdf5BottomSegment::hasParseUp() const {
        return getTopParseIndex() != NULL_INDEX;
    }

    inline hal_index_t Hdf5BottomSegment::getArrayIndex() const {
        return _index;
    }

    inline bool Hdf5BottomSegment::isFirst() const {
        assert(getSequence() != NULL);
        return _index == 0 || _index == (hal_index_t)getSequence()->getBottomSegmentArrayIndex();
    }

    inline bool Hdf5BottomSegment::isLast() const {
        assert(getSequence() != NULL);
        return _index == (hal_index_t)_array->getSize() - 1 ||
               _index == getSequence()->getBottomSegmentArrayIndex() + (hal_index_t)getSequence()->getNumBottomSegments() - 1;
    }

    inline hal_index_t Hdf5BottomSegment::getLeftChildIndex(hal_size_t i) const {
        assert(isFirst() == false);
        Hdf5BottomSegment leftSeg(getHdf5Genome(), _array, _index - 1);
        return leftSeg.getChildIndex(i);
    }

    inline hal_index_t Hdf5BottomSegment::getRightChildIndex(hal_size_t i) const {
        assert(isLast() == false);
        Hdf5BottomSegment rightSeg(dynamic_cast<Hdf5Genome *>(_genome), _array, _index + 1);
        return rightSeg.getChildIndex(i);
    }
}

#endif
// Local Variables:
// mode: c++
// End:
