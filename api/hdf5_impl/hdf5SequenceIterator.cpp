/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "hdf5SequenceIterator.h"
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;
using namespace H5;
using namespace hal;

Hdf5SequenceIterator::Hdf5SequenceIterator(Hdf5Genome *genome, hal_index_t index)
    : _sequence(genome, &genome->_sequenceIdxArray, &genome->_sequenceNameArray, index) {
}

Hdf5SequenceIterator::~Hdf5SequenceIterator() {
}

SequenceIteratorPtr Hdf5SequenceIterator::clone() const {
    Hdf5SequenceIterator *seqIt = new Hdf5SequenceIterator(_sequence._genome, _sequence._index);
    return SequenceIteratorPtr(seqIt);
}

void Hdf5SequenceIterator::toNext() {
    ++_sequence._index;
}

void Hdf5SequenceIterator::toPrev() {
    --_sequence._index;
}

bool Hdf5SequenceIterator::atEnd() const {
    return (_sequence._index < 0) or (_sequence._index >= (hal_index_t)_sequence._genome->_sequenceNameArray.getSize());
}

Sequence *Hdf5SequenceIterator::getSequence() {
    assert(_sequence._index >= 0 && _sequence._index < (hal_index_t)_sequence._genome->_sequenceNameArray.getSize());
    // don't return local sequence pointer.  give cached pointer from
    // genome instead (so it will not expire when iterator moves!)
    return _sequence._genome->getSequence(_sequence.getName());
}

bool Hdf5SequenceIterator::equals(SequenceIteratorPtr other) const {
    const Hdf5SequenceIterator *h5Other = reinterpret_cast<const Hdf5SequenceIterator *>(other.get());
    assert(_sequence.getGenome() == h5Other->_sequence.getGenome());
    return _sequence._index == h5Other->_sequence._index;
}
