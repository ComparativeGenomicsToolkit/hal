/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoTopSegment.h"
#include "halDnaIterator.h"
#include "yomoBottomSegment.h"
#include "yomoGenome.h"
#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;
using namespace H5;
using namespace hal;

const size_t YomoTopSegment::genomeIndexOffset = 0;
const size_t YomoTopSegment::bottomIndexOffset = sizeof(hal_index_t);
const size_t YomoTopSegment::parIndexOffset = bottomIndexOffset + sizeof(hal_index_t);
const size_t YomoTopSegment::parentIndexOffset = parIndexOffset + sizeof(hal_index_t);
const size_t YomoTopSegment::parentReversedOffset = parentIndexOffset + sizeof(hal_index_t);
const size_t YomoTopSegment::totalSize = parentReversedOffset + sizeof(bool);

YomoTopSegment::YomoTopSegment(YomoGenome *genome, YomoExternalArray *array, hal_index_t index)
    : TopSegment(genome, index), _array(array) {
}

void YomoTopSegment::setCoordinates(hal_index_t startPos, hal_size_t length) {
    if ((_genome != NULL) &&
        (startPos >= (hal_index_t)_genome->getSequenceLength() || startPos + length > _genome->getSequenceLength())) {
        throw hal_exception("Trying to set top segment coordinate out of range");
    }

    _array->setValue(_index, genomeIndexOffset, startPos);
    _array->setValue(_index + 1, genomeIndexOffset, startPos + length);
}

hal_offset_t YomoTopSegment::getBottomParseOffset() const {
    assert(_index >= 0);
    hal_offset_t offset = 0;
    hal_index_t bottomIndex = getBottomParseIndex();
    if (bottomIndex != NULL_INDEX) {
        YomoGenome *genome = dynamic_cast<YomoGenome *>(_genome);
        YomoBottomSegment bs(genome, &genome->_bottomArray, bottomIndex);
        assert(bs.getStartPosition() <= getStartPosition());
        assert((hal_index_t)(bs.getStartPosition() + bs.getLength()) >= getStartPosition());
        offset = getStartPosition() - bs.getStartPosition();
    }
    return offset;
}

bool YomoTopSegment::isCanonicalParalog() const {
    bool isCanon = false;
    if (hasParent()) {
        YomoGenome *parGenome = const_cast<YomoGenome *>(dynamic_cast<const YomoGenome *>(_genome->getParent()));

        YomoBottomSegment parent(parGenome, &parGenome->_bottomArray, getParentIndex());
        hal_index_t childGenomeIndex = parGenome->getChildIndex(_genome);
        isCanon = parent.getChildIndex(childGenomeIndex) == _index;
    }
    return isCanon;
}

void YomoTopSegment::print(std::ostream &os) const {
    os << "YOMO Top Segment";
}

// YOMO SPECIFIC
H5::CompType YomoTopSegment::dataType() {
    // the in-memory representations and yomo representations
    // don't necessarily have to be the same, but it simplifies
    // testing for now.
    assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
    assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
    assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
    assert(PredType::NATIVE_CHAR.getSize() == sizeof(bool));

    H5::CompType dataType(totalSize);
    dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
    dataType.insertMember("bottomIdx", bottomIndexOffset, PredType::NATIVE_INT64);
    dataType.insertMember("paralogyIdx", parIndexOffset, PredType::NATIVE_INT64);
    dataType.insertMember("parentIdx", parentIndexOffset, PredType::NATIVE_INT64);
    dataType.insertMember("reverseFlag", parentReversedOffset, PredType::NATIVE_CHAR);

    return dataType;
}
