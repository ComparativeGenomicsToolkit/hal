/*
 * Copyright (C) 2012 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2019 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef _HDF5SEQUENCEITERATOR_H
#define _HDF5SEQUENCEITERATOR_H

#include "halSequenceIterator.h"
#include "hdf5ExternalArray.h"
#include "hdf5Genome.h"
#include "hdf5Sequence.h"
#include <H5Cpp.h>

namespace hal {

    class Hdf5SequenceIterator : public SequenceIterator {
      public:
        Hdf5SequenceIterator(Hdf5Genome *genome, hal_index_t index);
        ~Hdf5SequenceIterator();

        // SEQUENCE ITERATOR METHODS
        SequenceIteratorPtr clone() const;
        void toNext();
        void toPrev();
        bool atEnd() const;
        const Sequence *getSequence() const;
        bool equals(SequenceIteratorPtr other) const;

      private:
        Hdf5Sequence _sequence;
    };
}
#endif
// Local Variables:
// mode: c++
// End:
