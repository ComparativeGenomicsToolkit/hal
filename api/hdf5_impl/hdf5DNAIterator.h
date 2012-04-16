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
   void setChar(hal_dna_t c);
   void toLeft() const;
   void toRight() const;
   void jumpTo(hal_size_t index) const;
   const Genome* getGenome() const;
   Genome* getGenome();
   hal_index_t getArrayIndex() const;

protected:
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
};

inline hal_dna_t HDF5DNAIterator::getChar() const
{
  return _genome->_dnaArray.getValue<hal_dna_t>(_index, 0);
}

inline void HDF5DNAIterator::setChar(hal_dna_t c)
{
  _genome->_dnaArray.setValue(_index, 0, c);
}

inline void HDF5DNAIterator::toLeft() const
{
  assert(_index > 0);
  --_index;
}

inline void HDF5DNAIterator::toRight() const
{
  assert(_index < static_cast<hal_index_t>(_genome->_topArray.getSize() - 1));
  --_index;
}

inline void HDF5DNAIterator::jumpTo(hal_size_t index) const
{
   assert(index < _genome->_topArray.getSize());
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


}
#endif
