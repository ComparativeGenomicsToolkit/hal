/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5DNACHARTYPE_H
#define _HDF5DNACHARTYPE_H

#include "rawH5ExternalArray.h"

namespace hal {

/**
 * Specialization of DNA character type for HDF5 external array. 
 */
struct DNACharType<RawH5ExternalArray>
{
   /** Return object with information for arrays containing top 
    * segments */
   static H5::PredType getTypeInfo();
};

// INLINE MEMBERS

inline H5::Datatype DNACharType<RawH5ExternalArray>::getTypeInfo()
{
  return PredType::NATIVE_CHAR;
}

}
#endif
