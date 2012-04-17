/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5DNAITERATOR_H
#define _HDF5DNAITERATOR_H

#include <cassert>
#include <H5Cpp.h>
#include "halDNAIterator.h"
#include "halCommon.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5DNA.h"

namespace hal {

class HDF5DNAIterator : public DNAIterator
{
public:
   
   HDF5DNAIterator(HDF5Genome* genome, hal_index_t index);
   ~HDF5DNAIterator();
   
   hal_dna_t getChar() const;
   hal_dna_t getCompChar() const;
   void setChar(hal_dna_t c);
   void toLeft() const;
   void toRight() const;
   void jumpTo(hal_size_t index) const;
   const Genome* getGenome() const;
   Genome* getGenome();
   hal_index_t getArrayIndex() const;

   void readString(std::string& outString, hal_size_t length,
                   hal_bool_t reversed = false) const;

   void writeString(const std::string& inString, hal_size_t length,
                    hal_bool_t reversed = false);

   inline bool inRange() const;

protected:
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};

inline bool HDF5DNAIterator::inRange() const
{
  return _index >= 0 && 
     _index < (hal_index_t)_genome->_dnaArray.getSize();
}

inline hal_dna_t HDF5DNAIterator::getChar() const
{
  assert(inRange() == true);
  return _genome->_dnaArray.getValue<hal_dna_t>(_index, 0);
}

inline hal_dna_t  HDF5DNAIterator::getCompChar() const
{
  assert(inRange() == true);
  return reverseComplement(getChar());
}

inline void HDF5DNAIterator::setChar(hal_dna_t c)
{
  if (inRange() == false) 
  {
    throw hal_exception("Trying to set character out of range");
  }
  _genome->_dnaArray.setValue(_index, 0, c);
}

inline void HDF5DNAIterator::toLeft() const
{
  --_index;
}

inline void HDF5DNAIterator::toRight() const
{
  ++_index;
}

inline void HDF5DNAIterator::jumpTo(hal_size_t index) const
{
  _index = static_cast<hal_index_t>(index);
}

inline const Genome* HDF5DNAIterator::getGenome() const
{
  return _genome;
}

inline Genome* HDF5DNAIterator::getGenome()
{
  return _genome;
}

inline hal_index_t HDF5DNAIterator::getArrayIndex() const
{
  return _index;
}

inline void HDF5DNAIterator::readString(std::string& outString,
                                        hal_size_t length, 
                                        hal_bool_t reversed) const
{
  assert(inRange() == true);
  outString.resize(length);
  if (reversed == false)
  {
    for (hal_size_t i = 0; i < length; ++i)
    {
      outString[i] = getChar();
      toRight();
    }
  }
  else
  {
    for (hal_index_t i = length - 1; i >= 0; --i)
    {
      outString[i] = getCompChar();
      toRight();
    }
  }
}

inline void HDF5DNAIterator::writeString(const std::string& inString,
                                         hal_size_t length,
                                         hal_bool_t reversed)
{
  assert(inRange() == true);
  if (reversed == false)
  {
    for (hal_size_t i = 0; i < length; ++i)
    {
      setChar(inString[i]);
      toRight();
    }
  }
  else
  {
    for (hal_index_t i = length - 1; i >= 0; --i)
    {
      setChar(reverseComplement(inString[i]));
      toRight();
    }
  }
}

}
#endif
