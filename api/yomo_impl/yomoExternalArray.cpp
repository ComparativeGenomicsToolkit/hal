/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "yomoExternalArray.h"
#include <cassert>
#include <iostream>

using namespace hal;
using namespace H5;
using namespace std;

/** Constructor */
YomoExternalArray::YomoExternalArray()
    : _file(NULL), _size(0), _chunkSize(0), _bufStart(0), _bufEnd(0), _bufSize(0), _buf(NULL), _dirty(false) {
}

/** Destructor */
YomoExternalArray::~YomoExternalArray() {
    delete[] _buf;
}

/* initialize the internal data buffer */
void YomoExternalArray::initBuf() {
    _bufSize = _chunkSize > 1 ? _chunkSize : _size;
    _bufStart = 0;
    _bufEnd = _bufSize - 1;
    delete[] _buf;
    _buf = new char[_bufSize * _dataSize];
}

// Create a new dataset in specifed location
void YomoExternalArray::create(PortableH5Location *file, const H5std_string &path, const DataType &dataType,
                               hsize_t numElements, const DSetCreatPropList *inCparms, hsize_t chunksInBuffer) {
    // copy in parameters
    _file = file;
    _path = path;
    _dataType = dataType;
    _size = numElements;
    _dataSize = _dataType.getSize();
    _dataSpace = DataSpace(1, &_size);

    DSetCreatPropList cparms;
    if (inCparms) {
        cparms.copy(*inCparms);
    }

    // resolve chunking size (0 = do not chunk)
    if (cparms.getLayout() == H5D_CHUNKED) {
        cparms.getChunk(1, &_chunkSize);
        // don't support chunking when size=1
        if (_size == 1) {
            cparms = DSetCreatPropList();
            _chunkSize = 0;
        }
        // something's wrong with the input chunk size.  do one chunk.
        else if (_chunkSize <= 1 || _chunkSize >= _size) {
            _chunkSize = _size;
            cparms.setChunk(1, &_chunkSize);
        }
        _chunkSize *= chunksInBuffer;
    } else {
        _chunkSize = 0;
    }

    // create the internal data buffer
    initBuf();

    // create the yomo array
    _dataSet = _file->createDataSet(_path, _dataType, _dataSpace, cparms);
    _chunkSpace = DataSpace(1, &_bufSize);
    assert(getSize() == numElements);
    assert(_bufSize > 0 || _size == 0);
}

// Load an existing dataset into memory
void YomoExternalArray::load(PortableH5Location *file, const H5std_string &path, hsize_t chunksInBuffer) {
    // load up the parameters
    _file = file;
    _path = path;
    _dataSet = _file->openDataSet(_path);
    _dataType = _dataSet.getDataType();
    _dataSpace = _dataSet.getSpace();
    _dataSize = _dataType.getSize();
    assert(_dataSpace.getSimpleExtentNdims() == 1);
    _dataSpace.getSimpleExtentDims(&_size, NULL);
    DSetCreatPropList cparms;
    cparms.copy(_dataSet.getCreatePlist());

    // resolve chunking size (0 = do not chunk)
    if (cparms.getLayout() == H5D_CHUNKED) {
        cparms.getChunk(1, &_chunkSize);
        _chunkSize *= chunksInBuffer;
        if (_chunkSize == 1) {
            throw hal_exception("YomoExternalArray::create: "
                                "chunkSize of 1 not supported");
        }
        if (_chunkSize > _size) {
            throw hal_exception("YomoExternalArray::create: "
                                "chunkSize > array size is not supported");
        }
    } else {
        _chunkSize = 0;
    }
    initBuf();
    _bufStart = _bufEnd + 1; // set out of range to ensure page happens
    _chunkSpace = DataSpace(1, &_bufSize);
    assert(_bufSize > 0 || _size == 0);
}

// Write the memory buffer back to the file
void YomoExternalArray::write() {
    if (_dirty) {
        _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
        _dataSet.write(_buf, _dataType, _chunkSpace, _dataSpace);
        _dirty = false;
    }
}

// Page chunk containing index i into memory
void YomoExternalArray::page(hsize_t i) {
    if (_dirty) {
        write();
    }
    _bufStart = (i / _bufSize) * _bufSize; // todo: review
    _bufEnd = _bufStart + _bufSize - 1;

    if (_bufEnd >= _size) {
        _bufEnd = _size - 1;
        _bufSize = _bufEnd - _bufStart + 1;
        _chunkSpace = DataSpace(1, &_bufSize);
    }

    _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
    _dataSet.read(_buf, _dataType, _chunkSpace, _dataSpace);
    _dirty = false;
    assert(_bufSize > 0 || _size == 0);
}
