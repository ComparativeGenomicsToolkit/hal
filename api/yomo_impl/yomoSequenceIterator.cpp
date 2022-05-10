/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoSequenceIterator.h"
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;
using namespace H5;
using namespace hal;

YomoSequenceIterator::YomoSequenceIterator(YomoGenome *genome, hal_index_t index)
    : _sequence(genome, &genome->_sequenceIdxArray, &genome->_sequenceNameArray, index) {
}

YomoSequenceIterator::~YomoSequenceIterator() {
}

SequenceIteratorPtr YomoSequenceIterator::clone() const {
    YomoSequenceIterator *seqIt = new YomoSequenceIterator(_sequence._genome, _sequence._index);
    return SequenceIteratorPtr(seqIt);
}

void YomoSequenceIterator::toNext() {
    ++_sequence._index;
}

void YomoSequenceIterator::toPrev() {
    --_sequence._index;
}

bool YomoSequenceIterator::atEnd() const {
    return (_sequence._index < 0) or (_sequence._index >= (hal_index_t)_sequence._genome->_sequenceNameArray.getSize());
}

Sequence *YomoSequenceIterator::getSequence() {
    assert(_sequence._index >= 0 && _sequence._index < (hal_index_t)_sequence._genome->_sequenceNameArray.getSize());
    // don't return local sequence pointer.  give cached pointer from
    // genome instead (so it will not expire when iterator moves!)
    return _sequence._genome->getSequence(_sequence.getName());
}

bool YomoSequenceIterator::equals(SequenceIteratorPtr other) const {
    const YomoSequenceIterator *h5Other = reinterpret_cast<const YomoSequenceIterator *>(other.get());
    assert(_sequence.getGenome() == h5Other->_sequence.getGenome());
    return _sequence._index == h5Other->_sequence._index;
}
