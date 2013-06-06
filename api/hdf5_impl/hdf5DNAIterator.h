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
   
   char getChar() const;
   void setChar(char c);
   void toLeft() const;
   void toRight() const;
   void jumpTo(hal_size_t index) const;
   void toReverse() const;
   bool getReversed() const;
   void setReversed(bool reversed) const;
   const Genome* getGenome() const;
   Genome* getGenome();
   const Sequence* getSequence() const;
   Sequence* getSequence();
   hal_index_t getArrayIndex() const;

   bool equals(DNAIteratorConstPtr& other) const;
   bool leftOf(DNAIteratorConstPtr& other) const;

   void readString(std::string& outString, hal_size_t length) const;

   void writeString(const std::string& inString, hal_size_t length);

   inline bool inRange() const;
   

protected:
   mutable hal_index_t _index;
   mutable HDF5Genome* _genome;
   mutable bool _reversed;
};

inline bool HDF5DNAIterator::inRange() const
{
  return _index >= 0 && 
     _index < (hal_index_t)_genome->_totalSequenceLength &&
     _index / 2 < (hal_index_t)_genome->_dnaArray.getSize();
}

inline char HDF5DNAIterator::getChar() const
{
  assert(inRange() == true);
  char c = _genome->_dnaArray.getValue<char>(_index / 2, 0);
  c = HDF5DNA::unpack(_index, c);
  if (_reversed)
  {
    c = reverseComplement(c);
  }
  return c;
}

inline void HDF5DNAIterator::setChar(char c)
{
  if (inRange() == false) 
  {
    throw hal_exception("Trying to set character out of range");
  }
  else if (isNucleotide(c) == false)
  {
    throw hal_exception(std::string("Trying to set invalid charachter: ") + c);
  }
  if (_reversed)
  {
    c = reverseComplement(c);
  }
  char* old = _genome->_dnaArray.getUpdate(_index / 2);
  HDF5DNA::pack(c, _index, (unsigned char&)*old);
  assert(getChar() == !_reversed ? c : reverseComplement(c));
}

inline void HDF5DNAIterator::toLeft() const
{
  _reversed ? ++_index : --_index;
}

inline void HDF5DNAIterator::toRight() const
{
  _reversed ? --_index : ++_index;
}

inline void HDF5DNAIterator::jumpTo(hal_size_t index) const
{
  _index = static_cast<hal_index_t>(index);
}

inline void HDF5DNAIterator::toReverse() const
{
  _reversed = !_reversed;
}

inline bool HDF5DNAIterator::getReversed() const
{
  return _reversed;
}

inline void HDF5DNAIterator::setReversed(bool reversed) const
{
  _reversed = reversed;
}

inline const Genome* HDF5DNAIterator::getGenome() const
{
  return _genome;
}

inline Genome* HDF5DNAIterator::getGenome()
{
  return _genome;
}

inline const Sequence* HDF5DNAIterator::getSequence() const
{
  return _genome->getSequenceBySite(_index);
}

inline Sequence* HDF5DNAIterator::getSequence()
{
  return _genome->getSequenceBySite(_index);
}

inline hal_index_t HDF5DNAIterator::getArrayIndex() const
{
  return _index;
}

inline bool HDF5DNAIterator::equals(DNAIteratorConstPtr& other) const
{
  const HDF5DNAIterator* h5Other = reinterpret_cast<
     const HDF5DNAIterator*>(other.get());
  assert(_genome == h5Other->_genome);
  return _index == h5Other->_index;
}

inline bool HDF5DNAIterator::leftOf(DNAIteratorConstPtr& other) const
{
  const HDF5DNAIterator* h5Other = reinterpret_cast<
     const HDF5DNAIterator*>(other.get());
  assert(_genome == h5Other->_genome);
  return _index < h5Other->_index;
}

inline void HDF5DNAIterator::readString(std::string& outString,
                                        hal_size_t length) const
{
  assert(length == 0 || inRange() == true);
  outString.resize(length);

  for (hal_size_t i = 0; i < length; ++i)
  {
    outString[i] = getChar();
    toRight();
  }
}

inline void HDF5DNAIterator::writeString(const std::string& inString,
                                         hal_size_t length)
{
  assert(length == 0 || inRange() == true);
  for (hal_size_t i = 0; i < length; ++i)
  {
    setChar(inString[i]);
    toRight();
  }
}

}
#endif
