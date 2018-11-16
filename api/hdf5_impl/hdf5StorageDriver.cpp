/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "hdf5StorageDriver.h"
#include "hdf5Genome.h"
#include "hdf5Alignment.h"

using namespace hal;

static const int UDC_FETCH_SIZE = 64 * 1024;  // size to bring in for UDC access

HDF5DnaAccess::HDF5DnaAccess(Hdf5Genome* genome,
                             Hdf5ExternalArray* dnaArray,
                             hal_index_t index):
    DnaAccess(0, 0, NULL),
    _genome(genome),
    _dnaArray(dnaArray) {
    fetch(index);
}

void HDF5DnaAccess::flush() {
    // ensure that marked dirty
    if (_dirty) {
        _dnaArray->setDirty();
    }
    _dirty = false;
}
    
void HDF5DnaAccess::fetch(hal_index_t index) const {
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
