/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _RAWH5EXTERNALARRAY_H
#define _RAWH5EXTERNALARRAY_H

#include <cassert>
#include <H5Cpp.h>
#include "halDefs.h"

namespace hal {

/** 
 * Wrapper for a 1-dimensional HDF5 array of fixed length.  Array objects
 * are defined (and typed) by the input datatype.  The array is paged into
 * memory chunk-by-chunk as needed (using the HDF5 cache as a back-end).  
 * We can't use compiler tpying of the input objects (and instead just 
 * expose the raw void* data) because the elements' sizes are not known
 * at compile time.  
 * Exposes same interface as RawH5Array (we skip virtual functions for 
 * performance reasons). Writing is implemented differently howerver,
 * where changes will be written on page-out rather than when write() is 
 * called.
 */
class RawH5ExternalArray
{
public:

   /** Constructor */
   RawH5ExternalArray();

   /** Destructor */
   virtual ~RawH5ExternalArray();

   /** Create a new dataset in specifed location*/
   void create(H5::H5File* file, const H5std_string& path,
              const H5::DataType& dataType, hsize_t numElements,
              hsize_t chunkSize);

   /** Load an existing dataset into memory*/
   void load(H5::H5File* file, const H5std_string& path);

   /** Write the memory buffer back to the file */
   void write();

   /** Access the raw data at given index */
   const void* get(hsize_t i);

   /** Write the raw data at given index */
   void* getUpdate(hsize_t i);
   
protected:

   /** Read chunk from file */
   void page(hsize_t i);

   /** Pointer to file that owns this dataset */
   H5::H5File* _file;
   /** Path of dataset in file */
   H5std_string _path;
   /** Datatype for array */
   H5::DataType _dataType;
   /** Dataspace (dimension information) of array */
   H5::DataSpace _dataSpace;
   /** The HDF5 array object */
   H5::DataSet _dataSet;
   /** Number of elements in the array (fixed length)*/
   hsize_t _numElements;
   /** Size of chunk in bytes */
   hsize_t _chunkSize;
   /** Size of datatype in bytes */
   hsize_t _dataSize;
   /** Index of first element in memory buffer */
   hsize_t _bufStart;
   /** Index of last element in memory buffer */
   hsize_t _bufEnd;
   /** Number of elements in memory buffer */
   hsize_t _bufLen;
   /** In-memory buffer */
   char* _buf;
   /** Dimensional information for in-memory buffer */
   H5::DataSpace _chunkSpace;
   /** Flag saying we should write to disk on write
    * or page-out calls (set by getUpdate()) */
   bool _dirty;

private:

   RawH5ExternalArray(const RawH5ExternalArray&);
   RawH5ExternalArray& operator=(const RawH5ExternalArray&);
};  

// INLINE MEMBERS

inline const void* RawH5ExternalArray::get(hsize_t i)
{
  assert(i < _numElements);
  if (i < _bufStart || i > _bufEnd)
  {
    page(i);
  }
  assert((i - _bufStart) * _dataSize < _chunkSize);
  return _buf + (i - _bufStart) * _dataSize;
}

inline void* RawH5ExternalArray::getUpdate(hsize_t i)
{
  assert(i < _numElements);
  if (i < _bufStart || i > _bufEnd)
  {
    page(i);
  }
  _dirty = true;
  assert((i - _bufStart) * _dataSize < _chunkSize);
  return _buf + (i - _bufStart) * _dataSize;
}

}
#endif
