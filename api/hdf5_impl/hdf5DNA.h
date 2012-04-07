/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5DNA_H
#define _HDF5DNA_H

#include <H5Cpp.h>
#include "halDefs.h"

namespace hal {

/** HDF5 datatype infromation for a DNA character
 */
class HDF5DNA
{
public:
   static H5::PredType dataType();
};

// inline members
inline H5::PredType HDF5DNA::dataType()
{
  return H5::PredType::NATIVE_CHAR;
}

}

#endif


