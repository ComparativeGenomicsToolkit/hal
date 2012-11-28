/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5DNA_H
#define _HDF5DNA_H

#include <cstdlib>
#include <H5Cpp.h>
#include "halDefs.h"

namespace hal {

/** HDF5 datatype infromation for a DNA character
 */
class HDF5DNA
{
public:
   static H5::PredType dataType();
   static char unpack(hal_index_t index, unsigned char packedChar);
   static void pack(char unpackedChar, hal_index_t index, 
                    unsigned char& packedChar);
};

// inline members
inline H5::PredType HDF5DNA::dataType()
{
  return H5::PredType::NATIVE_UINT8;
}

// we store two characters per byte. bit 1 is set for capital letter
// bits 2,3,4 determine character (acgtn).
inline char HDF5DNA::unpack(hal_index_t index, unsigned char packedChar)
{
  if (index % 2 == 0)
  {
    packedChar = packedChar >> 4;
  }
  bool capital = packedChar & 8U;
  char val;
  packedChar &= 7U;
  switch(packedChar)
  {
  case 0U : val = 'a'; break;
  case 1U : val = 'c'; break;
  case 2U : val = 'g'; break;
  case 3U : val = 't'; break;
  case 4U : val = 'n'; break;
  default : val = 'x'; break;
  }
  if (capital)
  {
    val = std::toupper(val);
  }
  return val;
}

inline void HDF5DNA::pack(char unpackedChar, hal_index_t index, 
                          unsigned char& packedChar)
{
  unsigned char val = 0U;
  if (std::isupper(unpackedChar))
  {
    val = 8U;
  }
  unpackedChar = std::tolower(unpackedChar);
  switch(unpackedChar) 
  {
  case 'a' : val |= 0U; break;
  case 'c' : val |= 1U; break;
  case 'g' : val |= 2U; break;
  case 't' : val |= 3U; break;
  case 'n' : val |= 4U; break;
  default : val |= 5U; break;
  }
  
  if (index % 2 == 0)
  {
    val = val << 4;
    packedChar &= 15U;
  }
  else
  {
    packedChar &= 240U;
  }
  packedChar |= val;
}

}

#endif


