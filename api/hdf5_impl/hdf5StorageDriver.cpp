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

HDF5DNAStorage::HDF5DNAStorage(HDF5Genome* genome,
                               HDF5ExternalArray* dnaArray,
                               hal_index_t index):
    DNAStorage(0, 0, NULL),
    _genome(genome),
    _dnaArray(dnaArray) {
    fetch(index);
}

void HDF5DNAStorage::fetch(hal_index_t index) const {
    if (_dirty) {
        _dnaArray->setDirty();
    }
    _dnaArray->page(index / 2);
    // DANGER _dnaArray is close-ended
    _startIndex = _dnaArray->getBufStart() / 2;
    _endIndex = (_dnaArray->getBufEnd() + 2) / 2; // 2 for open-ended convert and round up
    _buffer = _dnaArray->getBuf();
    _dirty = false;
}
