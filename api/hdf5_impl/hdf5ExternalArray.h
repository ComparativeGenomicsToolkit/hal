/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5EXTERNALARRAY_H
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
 * at compile time, and we don't want to move it around once its read.
 */
class HDF5ExternalArray
{
public:
 
   /** Constructor */
   HDF5ExternalArray();

   /** Destructor */
   virtual ~HDF5ExternalArray();

   /** Create a new dataset in specifed location
    * @param file Pointer to the HDF5 file in which to create array
    * @param path  Path within the file for the new array
    * @param dataType HDF5 Datatype describing contents of array
    * @param numElements Fixed length of the new array
    * @param cparms  HDF5 options (ie chunking and compression etc) 
    */
   void create(H5::H5File* file, 
               const H5std_string& path, 
               const H5::DataType& dataType,
               hsize_t numElements,
               const H5::DSetCreatPropList& cparms =
               H5::DSetCreatPropList::DEFAULT);
 
   /** Load an existing dataset into memory
     * @param file Pointer to the HDF5 file in which to create array
     * @param path  Path within the file for the new array
     * @param chunksInBuffer HDF5 options (ie chunking and compression etc) 
     * 0: load entire array into buffer
     * 1: use default chunking (from dataset)
     * N: buffersize will be N chunks. 
     */
   void load(H5::H5File* file, const H5std_string& path,
             hsize_t chunksInBuffer = 1);
   
   /** Write the memory buffer back to the file */
   void write();

   /** Access the raw data at given index
    * @param i index of element to retrieve for reading
    */
   const void* get(hsize_t i);

   /** Write the raw data at given index 
    * @param i index of element to retrieve for updating
    */
   void* getUpdate(hsize_t i);

   /** Number of elements in array */
   hsize_t getSize() const;
   
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
   hsize_t _size;
   /** Size of chunk in elements */
   hsize_t _chunkSize;
   /** Size of datatype in bytes */
   hsize_t _dataSize;
   /** Index of first element in memory buffer */
   hsize_t _bufStart;
   /** Index of last element in memory buffer */
   hsize_t _bufEnd;
   /** Number of elements in memory buffer */
   hsize_t _bufSize;
   /** In-memory buffer */
   char* _buf;
   /** Dimensional information for in-memory buffer */
   H5::DataSpace _chunkSpace;
   /** Flag saying we should write to disk on write
    * or page-out calls (set by getUpdate()) */
   bool _dirty;

private:

   HDF5ExternalArray(const HDF5ExternalArray&);
   HDF5ExternalArray& operator=(const HDF5ExternalArray&);
};  

// INLINE MEMBERS

inline const void* HDF5ExternalArray::get(hsize_t i)
{
  assert(i < _size);
  if (i < _bufStart || i > _bufEnd)
  {
    page(i);
  }
  assert((i - _bufStart) < _bufSize);
  return _buf + (i - _bufStart) * _dataSize;
}

inline void* HDF5ExternalArray::getUpdate(hsize_t i)
{
  assert(i < _size);
  if (i < _bufStart || i > _bufEnd)
  {
    page(i);
  }
  _dirty = true;
  assert((i - _bufStart) < _bufSize);
  return _buf + (i - _bufStart) * _dataSize;
}

inline hsize_t HDF5ExternalArray::getSize() const
{
  return _size;
}

}
#endif
