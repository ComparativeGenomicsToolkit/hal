/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <cassert>
#include <iostream>
#include "RawH5ExternalArray.h"

using namespace hal;
using namespace H5;

/** Constructor */
RawH5ExternalArray::RawH5ExternalArray() :
  _file(NULL),
  _numElements(0),
  _chunkSize(0),
  _dataSize(0),
  _bufStart(0),
  _bufEnd(0),
  _bufLen(0),
  _buf(0),
  _dirty(false)
{}

/** Destructor */
RawH5ExternalArray::~RawH5ExternalArray()
{
  delete [] _buf;
}

/** Create a new dataset in specifed location*/
void RawH5ExternalArray::create(H5File* file, 
                                const H5std_string& path,
                                const DataType& dataType,
                                hsize_t numElements,
                                hsize_t chunkSize)
{
  if (chunkSize % dataType.getSize() != 0)
  {
    throw (DataSetIException("RawH5ExternalArray::create",
                             "chunkSize not divisible by datatype size"));
  }
  assert(numElements > 0);
  
  _file = file;
  _path = path;
  _dataType = dataType;
  _numElements = numElements;
  _chunkSize = chunkSize;
  _dataSize = _dataType.getSize();
  _dataSpace = DataSpace(1, &_numElements);
  _dataSet = _file->createDataSet(_path, _dataType, _dataSpace);
  delete [] _buf;
  _buf = new char[_chunkSize];
  _bufLen = _chunkSize / _dataSize;  
  _bufStart = 0;
  _bufEnd = _bufStart + _bufLen - 1;
  _chunkSpace = DataSpace(1, &_bufLen);

  std::cout << "datasize = " << _dataSize
            << " chunksize = " << _chunkSize
            << " bufStart = " << _bufStart
            << " bufEnd = " << _bufEnd
            << " bufLen = " << _bufLen << std::endl;
}

/** Load an existing dataset into memory*/
void RawH5ExternalArray::load(H5File* file, 
                              const H5std_string& path)
{
  _file = file;
  _path = path;
  _dataSet = _file->openDataSet(_path);
  _dataType = _dataSet.getDataType();
  _dataSpace = _dataSet.getSpace();
  assert(_dataSpace.getSimpleExtentNdims() == 1);
  _dataSpace.getSimpleExtentDims(&_numElements, NULL);
  DSetCreatPropList cparms = _dataSet.getCreatePlist();
  assert(cparms.getLayout() == H5D_CHUNKED);
  cparms.getChunk(1, &_chunkSize);
  delete [] _buf;
  _buf = new char[_chunkSize];
  page(0);
}

/** Write the memory buffer back to the file */
void RawH5ExternalArray::write()
{
  if (_dirty == true)
  {
    _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufLen, &_bufStart);
    _dataSet.write(_buf, _dataType, _chunkSpace, _dataSpace);
  }
}

/** Page chunk containing index i into memory */
void RawH5ExternalArray::page(hsize_t i)
{
  if (_dirty == true)
  {
    write();
  }
  _bufStart = (i / _bufLen) * _bufLen; // todo: review
  _bufEnd = _bufStart + _bufLen - 1;
  _dataSpace.selectHyperslab(H5S_SELECT_SET, &_bufLen, &_bufStart);
  _dataSet.read(_buf, _dataType, _chunkSpace, _dataSpace);
  _dirty = false;
}
