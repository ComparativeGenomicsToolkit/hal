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
   static inline char unpack(hal_index_t index, unsigned char packedChar);
   static inline void pack(char unpackedChar, hal_index_t index, 
                           unsigned char& packedChar);
private:
    /* map of character to encoding for both upper and lower case */
    static constexpr uint8_t pack_map[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 4, 11, 4, 4, 4, 10, 4, 4, 4, 4, 4, 4, 12, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 3, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

    /* map of 8 bit encoding to character */
    static constexpr char unpack_map[16] = {'a', 't', 'g', 'c', 'n', '\x00', '\x00', '\x00', 'A', 'T', 'G', 'C', 'N', '\x00', '\x00', '\x00'};
};

// inline members
inline H5::PredType HDF5DNA::dataType()
{
  return H5::PredType::NATIVE_UINT8;
}

// we store two characters per byte (one per nibble). bit 1 is set for capital letter
// bits 2,3,4 determine character (acgtn).
inline char HDF5DNA::unpack(hal_index_t index, unsigned char packedChar)
{
    uint8_t code = (index & 1) ? (packedChar & 0x0F) : (packedChar >> 4);
    return unpack_map[code];
}

inline void HDF5DNA::pack(char unpackedChar, hal_index_t index, 
                          unsigned char& packedChar)
{
    uint8_t code = pack_map[uint8_t(unpackedChar)];
    packedChar = (index & 1) ? ((packedChar & 0xF0) | code) : ((packedChar & 0x0F) | (code << 4));
}

}

#endif


