/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include "hdf5ExternalArray.h"

using namespace hal;
using namespace H5;
using namespace std;

/** Constructor */
HDF5ExternalArray::HDF5ExternalArray() :
  _file(NULL),
  _size(0),
  _chunkSize(0),
  _bufStart(0),
  _bufEnd(0),
  _bufSize(0),
  _buf(NULL),
  _dirty(false)
{}

/** Destructor */
HDF5ExternalArray::~HDF5ExternalArray()
{
  delete [] _buf;
}

// Create a new dataset in specifed location
void HDF5ExternalArray::create(CommonFG* file, 
                               const H5std_string& path, 
                               const DataType& dataType,
                               hsize_t numElements,
                               const DSetCreatPropList* inCparms,
                               hsize_t chunksInBuffer)
{
  // copy in parameters
  _file = file;
  _path = path;
  _dataType = dataType;
  _size = numElements;
  _dataSize = _dataType.getSize();
  _dataSpace = DataSpace(1, &_size);
  
  DSetCreatPropList cparms;
  if (inCparms)
  {
    cparms.copy(*inCparms);
  }
  
  // resolve chunking size (0 = do not chunk)
  if (cparms.getLayout() == H5D_CHUNKED)
  {
    cparms.getChunk(1, &_chunkSize);
    // don't support chunking when size=1
    if (_size == 1)
    {
      cparms = DSetCreatPropList();
      _chunkSize = 0;
    }
    // something's wrong with the input chunk size.  do one chunk. 
    else if (_chunkSize <= 1 || _chunkSize >= _size)
    {
      _chunkSize = _size;
      cparms.setChunk(1, &_chunkSize);
    }
    _chunkSize *= chunksInBuffer;
  }
  else
  {
    _chunkSize = 0;
  }
  
  // create the internal data buffer
  _bufSize = _chunkSize > 1 ? _chunkSize : _size;  
  _bufStart = 0;
  _bufEnd = _bufStart + _bufSize - 1;
  delete [] _buf;
  _buf = new char[_bufSize * _dataSize];

  // create the hdf5 array
  _dataSet = _file->createDataSet(_path, _dataType, _dataSpace, cparms);
  _chunkSpace = DataSpace(1, &_bufSize);
  assert(getSize() == numElements);
  assert(_bufSize > 0 || _size == 0);
}

// Load an existing dataset into memory
void HDF5ExternalArray::load(CommonFG* file, const H5std_string& path,
                             hsize_t chunksInBuffer)
{
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
  if (cparms.getLayout() == H5D_CHUNKED)
  {
    cparms.getChunk(1, &_chunkSize);
    _chunkSize *= chunksInBuffer;
    if (_chunkSize == 1)
    {
      throw hal_exception("HDF5ExternalArray::create: " 
                          "chunkSize of 1 not supported");
    }
    if (_chunkSize > _size)
    {
      throw hal_exception("HDF5ExternalArray::create: "
                          "chunkSize > array size is not supported");
    }
  }
  else
  {
    _chunkSize = 0;
  }
  
  // create the internal data buffer
  _bufSize = _chunkSize > 1 ? _chunkSize : _size;  
  _bufEnd = _bufStart + _bufSize - 1;
  // set out of range to ensure page happens
  _bufStart = _bufEnd + 1;
  delete [] _buf;
  _buf = new char[_bufSize * _dataSize];

  assert(_bufSize > 0 || _size == 0);
}

// Write the memory buffer back to the file 
void HDF5ExternalArray::write()
{
  if (_dirty == true)
  {
    _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
    _dataSet.write(_buf, _dataType, _chunkSpace, _dataSpace);
  }
}

// Page chunk containing index i into memory 
void HDF5ExternalArray::page(hsize_t i)
{
  if (_dirty == true)
  {
    write();
  }
  _bufSize = _chunkSize > 1 ? _chunkSize : _size;  
  _bufStart = (i / _bufSize) * _bufSize; // todo: review
  _bufEnd = _bufStart + _bufSize - 1;  

  if (_bufEnd >= _size)
  {
    _bufEnd = _size - 1;
    _bufSize = _bufEnd - _bufStart + 1;
  }

  _chunkSpace = DataSpace(1, &_bufSize);
  _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufSize, &_bufStart);
  _dataSet.read(_buf, _dataType, _chunkSpace, _dataSpace);
  _dirty = false;
  assert(_bufSize > 0 || _size == 0);
}
