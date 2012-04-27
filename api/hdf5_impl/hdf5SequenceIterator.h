/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5SEQUENCEITERATOR_H
#define _HDF5SEQUENCEITERATOR_H

#include <H5Cpp.h>
#include "halSequenceIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5Sequence.h"

namespace hal {

class HDF5BottomSegmentIterator;

class HDF5SequenceIterator : public SequenceIterator
{
   friend class HDF5BottomSegmentIterator;
public:
   
   HDF5SequenceIterator(HDF5Genome* genome, hal_index_t index);
   ~HDF5SequenceIterator();
   
   // SEQUENCE ITERATOR METHODS
   SequenceIteratorPtr copy();
   SequenceIteratorConstPtr copy() const;
   void toNext() const;
   void toPrev() const;
   Sequence* getSequence();
   const Sequence* getSequence() const;
   bool equals(SequenceIteratorConstPtr other) const;

protected:
   HDF5Sequence _sequence;
};

}
#endif
