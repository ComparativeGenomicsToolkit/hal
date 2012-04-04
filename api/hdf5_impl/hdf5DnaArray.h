/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5DNAARRAY_H
#define _HDF5DNAARRAY_H

#include <cassert>
#include <H5Cpp.h>
#include "rawH5ExternalArray.h"

namespace hal {

/** 
 * Wraps the RawH5ExternalArray with interface tailored to storing and accessing
 * only DNA characters
 */
class Hdf5DnaArray
{
public:

   /** Constructor */
   Hdf5DnaArray();  

   /** Destructor */
   virtual ~Hdf5DnaArray();
   
   /** Create a new array (overloads method in parent)
    * @param file HDF5 file in which to add new array dataset
    * @param path location of new array in file
    * @param size Fixed length of the new array
    * @param cparams Creation parameters for new array (chunking, zipping) */
   void create(H5File* file,
               const std::string& path,
               hsize_t size, 
               const H5::DSetCreatPropList& cparms = 
               H5::DSetCreatPropList::DEFAULT);
   
   /** Open an existing array 
    * @param file HDF5 file containing array to open
    * @param path location of array in file */
   void open(H5File* file,
             const std::string& path);

   /** Write any unsaved buffer contents back to the file */
   void write();
             
   /** Get read/write iterator 
    * @param offset position of iterator in array */
   Hdf5DnaIterator getDnaIterator(hsize_t offset = 0);

   /** Get read-only iterator
    * @param offset position of iterator in array */
   Hdf5DnaConstIterator getDnaConstIterator(hsize_t offset = 0);

   /** Get size of array */
   hsize_t size();
   
protected:
   
   RawH5ExternalArray _array;
};

// INLINE METHODS

inline Hdf5DnaArray::getDnaIterator(hsize_t offset)
{
  assert(offset < size());
  return DnaIterator(_array, offset);
}

inline Hdf5DnaConstIterator(hsize_t offset)
{
  assert(offset < size());
  return DnaIterator(_array, offset);
}

}
