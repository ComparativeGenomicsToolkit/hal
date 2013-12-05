/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hdf5SequenceIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

HDF5SequenceIterator::HDF5SequenceIterator(HDF5Genome* genome, 
                                           hal_index_t index) :
_sequence(genome, &genome->_sequenceIdxArray, 
          &genome->_sequenceNameArray, index)
{
  
}

HDF5SequenceIterator::~HDF5SequenceIterator()
{

}
   
SequenceIteratorPtr HDF5SequenceIterator::copy()
{
  HDF5SequenceIterator* newIt = new HDF5SequenceIterator(
    _sequence._genome, _sequence._index);
  return SequenceIteratorPtr(newIt);
}

SequenceIteratorConstPtr HDF5SequenceIterator::copy() const
{
  HDF5SequenceIterator* newIt = new HDF5SequenceIterator(
    _sequence._genome, _sequence._index);
  return SequenceIteratorConstPtr(newIt);
}

void HDF5SequenceIterator:: toNext() const
{
  ++_sequence._index;
}

void HDF5SequenceIterator::toPrev() const
{
  --_sequence._index;
}

Sequence* HDF5SequenceIterator::getSequence()
{
  assert(_sequence._index >= 0 && _sequence._index < 
         (hal_index_t)_sequence._genome->_sequenceNameArray.getSize()); 
  // don't return local sequence pointer.  give cached pointer from
  // genome instead (so it will not expire when iterator moves!)
  return _sequence._genome->getSequence(_sequence.getName());
}

const Sequence* HDF5SequenceIterator::getSequence() const
{
  assert(_sequence._index >= 0 && _sequence._index < 
         (hal_index_t)_sequence._genome->_sequenceNameArray.getSize());
  // don't return local sequence pointer.  give cached pointer from
  // genome instead (so it will not expire when iterator moves!)
  return _sequence._genome->getSequence(_sequence.getName());
}

bool HDF5SequenceIterator::equals(SequenceIteratorConstPtr other) const
{
  const HDF5SequenceIterator* h5Other = reinterpret_cast<
     const HDF5SequenceIterator*>(other.get());
  assert(_sequence.getGenome() == h5Other->_sequence.getGenome());
  return _sequence._index == h5Other->_sequence._index;
}

