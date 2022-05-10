/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoSequence.h"
#include "halBottomSegmentIterator.h"
#include "halColumnIterator.h"
#include "halDnaIterator.h"
#include "halGappedBottomSegmentIterator.h"
#include "halGappedTopSegmentIterator.h"
#include "halRearrangement.h"
#include "halTopSegmentIterator.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <set>
#include <string>

using namespace std;
using namespace H5;
using namespace hal;

const size_t YomoSequence::startOffset = 0;
const size_t YomoSequence::topSegmentArrayIndexOffset = sizeof(hal_size_t);
const size_t YomoSequence::bottomSegmentArrayIndexOffset = topSegmentArrayIndexOffset + sizeof(hal_size_t);
const size_t YomoSequence::totalSize = bottomSegmentArrayIndexOffset + sizeof(hal_size_t);

YomoSequence::YomoSequence(YomoGenome *genome, YomoExternalArray *idxArray, YomoExternalArray *nameArray, hal_index_t index)
    : _idxArray(idxArray), _nameArray(nameArray), _index(index), _genome(genome) {
}

YomoSequence::~YomoSequence() {
}

H5::CompType YomoSequence::idxDataType() {
    // the in-memory representations and yomo representations
    // don't necessarily have to be the same, but it simplifies
    // testing for now.
    assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
    CompType dataType(totalSize);
    dataType.insertMember("start", startOffset, PredType::NATIVE_HSIZE);
    dataType.insertMember("topSegmentArrayIndexOffset", topSegmentArrayIndexOffset, PredType::NATIVE_HSIZE);
    dataType.insertMember("bottomSegmentArrayIndexOffset", bottomSegmentArrayIndexOffset, PredType::NATIVE_HSIZE);
    return dataType;
}

H5::StrType YomoSequence::nameDataType(hal_size_t maxNameLength) {
    assert(PredType::NATIVE_HSIZE.getSize() == sizeof(hal_size_t));
    assert(PredType::NATIVE_CHAR.getSize() == sizeof(char));
    StrType strType(PredType::NATIVE_CHAR, (maxNameLength + 1) * sizeof(char));
    return strType;
}

// SEQUENCE INTERFACE
string YomoSequence::getName() const {
    return _nameArray->get(_index);
}

string YomoSequence::getFullName() const {
    assert(_genome != NULL);
    return _genome->getName() + '.' + getName();
}

const Genome *YomoSequence::getGenome() const {
    return _genome;
}

Genome *YomoSequence::getGenome() {
    return _genome;
}

hal_index_t YomoSequence::getStartPosition() const {
    return _idxArray->getValue<hal_size_t>(_index, startOffset);
}

hal_index_t YomoSequence::getEndPosition() const {
    return _idxArray->getValue<hal_size_t>(_index + 1, startOffset) - 1;
}

hal_index_t YomoSequence::getArrayIndex() const {
    return _index;
}

hal_index_t YomoSequence::getTopSegmentArrayIndex() const {
    return (hal_index_t)_idxArray->getValue<hal_size_t>(_index, topSegmentArrayIndexOffset);
}

hal_index_t YomoSequence::getBottomSegmentArrayIndex() const {
    return (hal_index_t)_idxArray->getValue<hal_size_t>(_index, bottomSegmentArrayIndexOffset);
}

// SEGMENTED SEQUENCE INTERFACE

hal_size_t YomoSequence::getSequenceLength() const {
    hal_size_t len = _idxArray->getValue<hal_size_t>(_index, startOffset);
    hal_size_t nlen = _idxArray->getValue<hal_size_t>(_index + 1, startOffset);
    assert(nlen >= len);
    return (nlen - len);
}

hal_size_t YomoSequence::getNumTopSegments() const {
    hal_index_t idx = _idxArray->getValue<hal_size_t>(_index, topSegmentArrayIndexOffset);
    hal_index_t nextIdx = _idxArray->getValue<hal_size_t>(_index + 1, topSegmentArrayIndexOffset);
    assert(nextIdx >= idx);
    return nextIdx - idx;
}

hal_size_t YomoSequence::getNumBottomSegments() const {
    hal_index_t idx = _idxArray->getValue<hal_size_t>(_index, bottomSegmentArrayIndexOffset);
    hal_index_t nextIdx = _idxArray->getValue<hal_size_t>(_index + 1, bottomSegmentArrayIndexOffset);
    assert(nextIdx >= idx);
    return nextIdx - idx;
}

TopSegmentIteratorPtr YomoSequence::getTopSegmentIterator(hal_index_t position) {
    hal_size_t idx = position + getTopSegmentArrayIndex();
    return _genome->getTopSegmentIterator(idx);
}

TopSegmentIteratorPtr YomoSequence::getTopSegmentIterator(hal_index_t position) const {
    hal_size_t idx = position + getTopSegmentArrayIndex();
    return _genome->getTopSegmentIterator(idx);
}

BottomSegmentIteratorPtr YomoSequence::getBottomSegmentIterator(hal_index_t position) {
    hal_size_t idx = position + getBottomSegmentArrayIndex();
    return _genome->getBottomSegmentIterator(idx);
}

BottomSegmentIteratorPtr YomoSequence::getBottomSegmentIterator(hal_index_t position) const {
    hal_size_t idx = position + getBottomSegmentArrayIndex();
    return _genome->getBottomSegmentIterator(idx);
}

DnaIteratorPtr YomoSequence::getDnaIterator(hal_index_t position) {
    return _genome->getDnaIterator(position + getStartPosition());
}

DnaIteratorPtr YomoSequence::getDnaIterator(hal_index_t position) const {
    return const_cast<YomoSequence *>(this)->getDnaIterator(position);
}

ColumnIteratorPtr YomoSequence::getColumnIterator(const std::set<const Genome *> *targets, hal_size_t maxInsertLength,
                                                  hal_index_t position, hal_index_t lastPosition, bool noDupes,
                                                  bool noAncestors, bool reverseStrand, bool unique, bool onlyOrthologs) const {
    hal_index_t idx = (hal_index_t)(position + getStartPosition());
    hal_index_t lastIdx;
    if (lastPosition == NULL_INDEX) {
        lastIdx = (hal_index_t)(getStartPosition() + getSequenceLength() - 1);
    } else {
        lastIdx = (hal_index_t)(lastPosition + getStartPosition());
    }
    if (position < 0 || lastPosition >= (hal_index_t)(getStartPosition() + getSequenceLength())) {
        throw hal_exception("YomoSequence::getColumnIterators: input indices (" + std::to_string(position) + ", " +
                            std::to_string(lastPosition) + ") out of bounds");
    }
    ColumnIterator *colIt = new ColumnIterator(getGenome(), targets, idx, lastIdx, maxInsertLength, noDupes, noAncestors,
                                               reverseStrand, unique, onlyOrthologs);
    return ColumnIteratorPtr(colIt);
}

void YomoSequence::getString(std::string &outString) const {
    getSubString(outString, 0, getSequenceLength());
}

void YomoSequence::setString(const std::string &inString) {
    setSubString(inString, 0, getSequenceLength());
}

void YomoSequence::getSubString(std::string &outString, hal_size_t start, hal_size_t length) const {
    outString.resize(length);
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->readString(outString, length);
}

void YomoSequence::setSubString(const std::string &inString, hal_size_t start, hal_size_t length) {
    if (length != inString.length()) {
        throw hal_exception("setString: input string of length " + std::to_string(inString.length()) +
                            " has length different from target string in sequence " + getName() + " which is of length " +
                            std::to_string(length));
    }
    DnaIteratorPtr dnaIt(getDnaIterator(start));
    dnaIt->writeString(inString, length);
}

RearrangementPtr YomoSequence::getRearrangement(hal_index_t position, hal_size_t gapLengthThreshold, double nThreshold,
                                                bool atomic) const {
    TopSegmentIteratorPtr top = getTopSegmentIterator(position);
    Rearrangement *rea = new Rearrangement(getGenome(), gapLengthThreshold, nThreshold, atomic);
    rea->identifyFromLeftBreakpoint(top);
    return RearrangementPtr(rea);
}

GappedTopSegmentIteratorPtr YomoSequence::getGappedTopSegmentIterator(hal_index_t i, hal_size_t gapThreshold,
                                                                      bool atomic) const {
    TopSegmentIteratorPtr topSegIt = getTopSegmentIterator(i);
    GappedTopSegmentIterator *gapTopSegIt = new GappedTopSegmentIterator(topSegIt, gapThreshold, atomic);
    return GappedTopSegmentIteratorPtr(gapTopSegIt);
}

GappedBottomSegmentIteratorPtr YomoSequence::getGappedBottomSegmentIterator(hal_index_t i, hal_size_t childIdx,
                                                                            hal_size_t gapThreshold, bool atomic) const {
    BottomSegmentIteratorPtr botSegIt = getBottomSegmentIterator(i);
    GappedBottomSegmentIterator *gapBotSegIt = new GappedBottomSegmentIterator(botSegIt, childIdx, gapThreshold, atomic);
    return GappedBottomSegmentIteratorPtr(gapBotSegIt);
}
// LOCAL

void YomoSequence::set(hal_size_t startPosition, const Sequence::Info &sequenceInfo, hal_size_t topSegmentStartIndex,
                       hal_size_t bottomSegmentStartIndex) {
    _idxArray->setValue(_index, startOffset, startPosition);
    _idxArray->setValue(_index + 1, startOffset, startPosition + sequenceInfo._length);
    _idxArray->setValue(_index, topSegmentArrayIndexOffset, topSegmentStartIndex);
    _idxArray->setValue(_index, bottomSegmentArrayIndexOffset, bottomSegmentStartIndex);
    _idxArray->setValue(_index + 1, topSegmentArrayIndexOffset, topSegmentStartIndex + sequenceInfo._numTopSegments);
    _idxArray->setValue(_index + 1, bottomSegmentArrayIndexOffset, bottomSegmentStartIndex + sequenceInfo._numBottomSegments);
    char *arrayBuffer = _nameArray->getUpdate(_index);
    strcpy(arrayBuffer, sequenceInfo._name.c_str());

    assert(getStartPosition() == (hal_index_t)startPosition);
    assert(getNumTopSegments() == sequenceInfo._numTopSegments);
    assert(getNumBottomSegments() == sequenceInfo._numBottomSegments);
    assert(getSequenceLength() == sequenceInfo._length);
}

// These functions look dangerous.  Dont think they're used.
void YomoSequence::setNumTopSegments(hal_size_t numTopSegments) {
    _idxArray->setValue(_index + 1, topSegmentArrayIndexOffset, getTopSegmentArrayIndex() + numTopSegments);
}

void YomoSequence::setNumBottomSegments(hal_size_t numBottomSegments) {
    _idxArray->setValue(_index + 1, bottomSegmentArrayIndexOffset, getBottomSegmentArrayIndex() + numBottomSegments);
}

void YomoSequence::setTopSegmentArrayIndex(hal_size_t topIndex) {
    _idxArray->setValue(_index, topSegmentArrayIndexOffset, topIndex);
}

void YomoSequence::setBottomSegmentArrayIndex(hal_size_t bottomIndex) {
    _idxArray->setValue(_index, bottomSegmentArrayIndexOffset, bottomIndex);
}

void YomoSequence::setName(const string &newName) {
    _genome->renameSequence(getName(), _index, newName);
}
