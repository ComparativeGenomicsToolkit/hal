/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "hdf5BottomSegment.h"
#include "halDnaIterator.h"
#include "hdf5Genome.h"
#include "hdf5TopSegment.h"
#include <cstdlib>
#include <string>

using namespace std;
using namespace H5;
using namespace hal;

const size_t Hdf5BottomSegment::genomeIndexOffset = 0;
const size_t Hdf5BottomSegment::lengthOffset = sizeof(hal_index_t);
const size_t Hdf5BottomSegment::topIndexOffset = lengthOffset + sizeof(hal_size_t);
const size_t Hdf5BottomSegment::firstChildOffset = topIndexOffset + sizeof(hal_index_t);
const size_t Hdf5BottomSegment::totalSize(hal_size_t numChildren) {
    return firstChildOffset + numChildren * (sizeof(hal_index_t) + sizeof(bool));
}

Hdf5BottomSegment::Hdf5BottomSegment(Hdf5Genome *genome, Hdf5ExternalArray *array, hal_index_t index)
    : BottomSegment(genome, index), _array(array) {
}

hal_size_t Hdf5BottomSegment::numChildrenFromDataType(const H5::DataType &dataType) {
    return (dataType.getSize() - firstChildOffset) / (sizeof(hal_index_t) + sizeof(bool));
}

hal_offset_t Hdf5BottomSegment::getTopParseOffset() const {
    assert(_index >= 0);
    hal_offset_t offset = 0;
    hal_index_t topIndex = getTopParseIndex();
    if (topIndex != NULL_INDEX) {
        Hdf5Genome *genome = dynamic_cast<Hdf5Genome *>(_genome);
        Hdf5TopSegment ts(genome, &genome->_topArray, topIndex);
        assert(ts.getStartPosition() <= getStartPosition());
        assert((hal_index_t)(ts.getStartPosition() + ts.getLength()) >= getStartPosition());
        offset = getStartPosition() - ts.getStartPosition();
    }
    return offset;
}

void Hdf5BottomSegment::setCoordinates(hal_index_t startPos, hal_size_t length) {
    assert(_index >= 0);
    if (_genome &&
        (startPos >= (hal_index_t)_genome->getSequenceLength() || startPos + length > _genome->getSequenceLength())) {
        throw hal_exception("Trying to set bottom segment coordinate out of range");
    }

    _array->setValue((hsize_t)_index, genomeIndexOffset, startPos);
    _array->setValue(_index + 1, genomeIndexOffset, startPos + length);
}

void Hdf5BottomSegment::print(std::ostream &os) const {
    os << "HDF5 Bottom Segment";
}

// HDF5 SPECIFIC
H5::CompType Hdf5BottomSegment::dataType(hal_size_t numChildren) {
    // the in-memory representations and hdf5 representations
    // don't necessarily have to be the same, but it simplifies
    // testing for now.
    assert(PredType::NATIVE_INT64.getSize() == sizeof(hal_index_t));
    assert(PredType::NATIVE_UINT64.getSize() == sizeof(hal_offset_t));
    assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
    assert(PredType::NATIVE_CHAR.getSize() == sizeof(bool));

    H5::CompType dataType(totalSize(numChildren));
    dataType.insertMember("genomeIdx", genomeIndexOffset, PredType::NATIVE_INT64);
    dataType.insertMember("length", lengthOffset, PredType::NATIVE_HSIZE);
    dataType.insertMember("topIdx", topIndexOffset, PredType::NATIVE_INT64);

    // note: I'm refactoring the old explicit representation into an array
    // with the same layout interleaving index and reverse flag. this
    // preserves compatibility with the original layout that put
    // everything as a field in "dataType" but exceeded hdf5's header
    // limit when numChildren > 500
    // https://github.com/ComparativeGenomicsToolkit/cactus/issues/477
    // https://github.com/ComparativeGenomicsToolkit/hal/issues/212
    if (numChildren > 0) {
        H5::CompType childIdxType(PredType::NATIVE_INT64.getSize() + PredType::NATIVE_CHAR.getSize());
        childIdxType.insertMember("childIdx", 0, PredType::NATIVE_INT64);
        childIdxType.insertMember("reverseFlag", PredType::NATIVE_INT64.getSize(), PredType::NATIVE_CHAR);
        hsize_t dims = numChildren;
        ArrayType childIdxArrayType(childIdxType, 1, &dims);
        dataType.insertMember("childIndexes", firstChildOffset, childIdxArrayType);
    }
    
    return dataType;
}
