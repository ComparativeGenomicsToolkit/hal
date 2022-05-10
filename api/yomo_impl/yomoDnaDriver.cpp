/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "yomoDnaDriver.h"
#include "yomoAlignment.h"
#include "yomoGenome.h"

using namespace hal;

YOMODnaAccess::YOMODnaAccess(YomoGenome *genome, YomoExternalArray *dnaArray, hal_index_t index)
    : DnaAccess(0, 0, NULL), _dnaArray(dnaArray) {
    fetch(index);
}

void YOMODnaAccess::flush() {
    // ensure that marked dirty
    if (_dirty) {
        _dnaArray->setDirty();
    }
    _dirty = false;
}

void YOMODnaAccess::fetch(hal_index_t index) const {
    if (_dirty) {
        _dnaArray->setDirty();
    }
    _dnaArray->page(index / 2);
    // DANGER _dnaArray is close-ended
    _startIndex = 2 * _dnaArray->getBufStart();
    _endIndex = 2 * (_dnaArray->getBufEnd() + 1);
    _buffer = _dnaArray->getBuf();
    _dirty = false;
}
