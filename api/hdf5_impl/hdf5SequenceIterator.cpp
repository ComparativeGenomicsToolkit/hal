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
#include "hdf5BottomSegmentIterator.h"

using namespace std;
using namespace H5;
using namespace hal;

HDF5SequenceIterator::HDF5SequenceIterator(HDF5Genome* genome, 
                                           hal_index_t index) :
  _sequence(genome, &genome->_sequenceArray, index)
{
}

HDF5SequenceIterator::~HDF5SequenceIterator()
{
}
   
SequenceIteratorPtr HDF5SequenceIterator::copy()
{
  return SequenceIteratorPtr();
}

SequenceIteratorConstPtr HDF5SequenceIterator::copy() const
{
  return SequenceIteratorConstPtr();
}

void HDF5SequenceIterator:: toNext() const
{

}

void HDF5SequenceIterator::toPrev() const
{

}

Sequence* HDF5SequenceIterator::getSequence()
{
  return NULL;
}

const Sequence* HDF5SequenceIterator::getSequence() const
{
  return NULL;
}

