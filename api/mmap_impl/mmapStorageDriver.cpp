/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "mmapStorageDriver.h"
#include "mmapGenome.h"
#include "mmapAlignment.h"

using namespace hal;

static const int UDC_FETCH_SIZE = 64 * 1024;  // size to bring in for UDC access

MMapDNAStorage::MMapDNAStorage(MMapGenome* genome,
                               hal_index_t index):
    DNAStorage(0, 0, NULL),
    _genome(genome),
    _isUdcProtocol(dynamic_cast<MMapAlignment*>(_genome->getAlignment())->getMMapFile()->isUdcProtocol()) {
    if (_isUdcProtocol) {
        fetch(index);
    } else {
        // for local mmap, just include the whole thing
        _endIndex = _genome->getTotalSequenceLength();
        _buffer = _genome->getDNA(0, 1);
    }
}


void MMapDNAStorage::fetch(hal_index_t index) const {
    if (_isUdcProtocol) {
        _startIndex = 2 * ( index / 2);  // even boundary
        _endIndex = std::max(hal_size_t(_startIndex + UDC_FETCH_SIZE),
                             _genome->getTotalSequenceLength());
        _buffer = _genome->getDNA(_startIndex / 2, (((_endIndex - _startIndex) + 1) / 2));
    } else {
        assert(false);  // this should never be called for local
    }
    _dirty = false;  // keep consistent, but not actually used
}
