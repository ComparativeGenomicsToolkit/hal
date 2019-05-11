/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "halMaskExtractor.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace hal;

MaskExtractor::MaskExtractor() {
}

MaskExtractor::~MaskExtractor() {
}

void MaskExtractor::extract(const Alignment *alignment, const Genome *genome, ostream *bedStream, hal_size_t extend,
                            double extendPct) {
    _alignment = AlignmentConstPtr(alignment);
    _genome = genome;
    _bedStream = bedStream;
    _extend = extend;
    _extendPct = extendPct;

    for (SequenceIteratorPtr seqIt = _genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        _sequence = seqIt->getSequence();
        if (_sequence->getSequenceLength() > 0) {
            _posCache.clear();
            addMaskedBasesToCache();
            extendCachedIntervals();
            writeCachedIntervals();
        }
    }
}

void MaskExtractor::addMaskedBasesToCache() {
    assert(_posCache.size() == 0);
    for (DnaIteratorPtr dna = _sequence->getDnaIterator(); !dna->atEnd(); dna->toRight()) {
        if (isMasked(dna->getBase())) {
            _posCache.insert(dna->getArrayIndex());
        }
    }
}

void MaskExtractor::extendCachedIntervals() {
    if (_extend == 0 && _extendPct == 0.) {
        return;
    }
    assert(_extend == 0 || _extendPct == 0.);

    hal_size_t start = (hal_size_t)_sequence->getStartPosition();
    hal_size_t last = (hal_size_t)(start + _sequence->getSequenceLength()) - 1;

    const PositionCache::IntervalSet *intervalSet = _posCache.getIntervalSet();
    PositionCache::IntervalSet::const_iterator i;
    vector<hal_index_t> padding;
    for (i = intervalSet->begin(); i != intervalSet->end(); ++i) {
        hal_size_t len = (hal_size_t)(i->first - i->second) + 1;
        hal_size_t pad = _extend ? _extend : (hal_size_t)(_extendPct * len);
        hal_size_t newFirst = max(start, i->second - pad);
        for (hal_size_t j = newFirst; j < (hal_size_t)i->second; ++j) {
            padding.push_back(j);
        }
        hal_size_t newLast = min(last, i->first + pad);
        for (hal_size_t j = i->first + 1; j <= newLast; ++j) {
            padding.push_back(j);
        }
    }
    for (hal_size_t k = 0; k < padding.size(); ++k) {
        _posCache.insert(padding[k]);
    }
}

void MaskExtractor::writeCachedIntervals() {
    hal_index_t start = _sequence->getStartPosition();
    const PositionCache::IntervalSet *intervalSet = _posCache.getIntervalSet();
    PositionCache::IntervalSet::const_iterator i;
    for (i = intervalSet->begin(); i != intervalSet->end(); ++i) {
        *_bedStream << _sequence->getName() << '\t' << i->second - start << '\t' << (i->first + 1) - start << '\n';
    }
}
